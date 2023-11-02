def estimate_output( wildcards ):
    output = list()
    for sample in config["SAMPLES"]:
        output.append( f"results/assembly/{sample}.annotated.gff" )
    if config["QC"]:
        output.append( "results/reports/qc_report.html" )
    return output


def format_cl_arg( wildcards, parameter, name, quoted=False ):
    if parameter:
        if quoted:
            return f"{name} {parameter:q}"
        else:
            return f"{name} {parameter}"
    else:
        return ""


rule all:
    input:
        estimate_output


rule quality_trimming:
    message:
        """Trim low quality bases from raw sequencing reads of {wildcards.sample}, using a window size of {params.window_size}, a minimum 
        average quality score of {params.window_minimum_quality}, and a minimum length of {params.minimum_length}.
        """
    input:
        reads1=lambda wildcards: config["SAMPLES"][wildcards.sample]["read1"],
        reads2=lambda wildcards: config["SAMPLES"][wildcards.sample]["read2"],
    params:
        window_size=config["quality_trimming"]["window_size"],
        window_minimum_quality=config["quality_trimming"]["minimum_quality"],
        minimum_length=config["quality_trimming"]["minimum_length"],
        baseout="intermediates/illumina/trimmed/{sample}.fastq.gz"
    output:
        read1_trimmed="intermediates/illumina/trimmed/{sample}_1P.fastq.gz",
        read1_unpaired=temp( "intermediates/illumina/trimmed/{sample}_1U.fastq.gz" ),
        read2_trimmed="intermediates/illumina/trimmed/{sample}_2P.fastq.gz",
        read2_unpaired=temp( "intermediates/illumina/trimmed/{sample}_2U.fastq.gz" ),
        stats="intermediates/illumina/trimmed/{sample}_stats.txt"
    threads: 4
    shell:
        """
        trimmomatic PE \
            -threads {threads} \
            {input.reads1} {input.reads2} \
            -baseout {params.baseout} \
            SLIDINGWINDOW:{params.window_size}:{params.window_minimum_quality} \
            MINLEN:{params.minimum_length} &> {output.stats}
        """


# Might be superflouous
rule adapter_trimming:
    message: "Trim Illumina sequencing adapters and PhiX control from the reads for {wildcards.sample}"
    input:
        read1=rules.quality_trimming.output.read1_trimmed,
        read2=rules.quality_trimming.output.read2_trimmed
    params:
        firstpass_reference=config["adapter_trimming"]["firstpass_reference"],
        firstpass_args=config["adapter_trimming"]["firstpass_arguments"],
        secondpass_reference=config["adapter_trimming"]["secondpass_reference"],
        secondpass_args=config["adapter_trimming"]["secondpass_arguments"]
    output:
        repaired1=temp( "intermediates/illumina/adapters_removed/{sample}_R1.repaired.fastq.gz" ),
        repaired2=temp( "intermediates/illumina/adapters_removed/{sample}_R2.repaired.fastq.gz" ),
        firstpass_read1=temp( "intermediates/illumina/adapters_removed/{sample}_R1.firstpass.fastq.gz" ),
        firstpass_read2=temp( "intermediates/illumina/adapters_removed/{sample}_R2.firstpass.fastq.gz" ),
        firstpass_log="intermediates/illumina/adapters_removed/{sample}.firstpass.log",
        cleaned_read1="intermediates/illumina/adapters_removed/{sample}_R1.cleaned.fastq.gz",
        cleaned_read2="intermediates/illumina/adapters_removed/{sample}_R2.cleaned.fastq.gz",
        secondpass_log="intermediates/illumina/adapters_removed/{sample}.cleaned.log"
    threads: 8
    shell:
        """
        repair.sh \
            in1={input.read1} in2={input.read2} \
            out1={output.repaired1} out2={output.repaired2} &&\
        bbduk.sh \
            in1={output.repaired1} in2={output.repaired2} \
            out1={output.firstpass_read1} out2={output.firstpass_read2} \
            ref={params.firstpass_reference} \
            stats={output.firstpass_log} \
            {params.firstpass_args} &&\
        bbduk.sh \
            in1={output.firstpass_read1} in2={output.firstpass_read2} \
            out1={output.cleaned_read1} out2={output.cleaned_read2} \
            ref={params.secondpass_reference} \
            stats={output.secondpass_log} \
            {params.secondpass_args}
        """


rule denovo_assembly:
    message: "Assemble reads using a wrapper around SPAdes for {wildcards.sample}"
    input:
        read1=rules.adapter_trimming.output.cleaned_read1,
        read2=rules.adapter_trimming.output.cleaned_read2
    params:
        assembler=config["assembly"]["assembler"],
        assembler_arguments=lambda wildcards: format_cl_arg(
            wildcards,config["assembly"]["assembler_arguments"],"--opts"
        ),
        minimum_contig_length=config["assembly"]["minimum_contig_length"],
        trim="--trim" if config["assembly"]["trim"] else "",
        disable_correction="" if config["assembly"]["correction"] else "--nocorr",
        disable_read_correction="" if config["assembly"]["read_correction"] else "--noreadcorr",
        disable_read_stitching="" if config["assembly"]["read_stitching"] else "--nostitch",
        temp_contigs=temp( "intermediates/illumina/assembly/{sample}/contigs.fa" ),
        temp_graph=temp( "intermediates/illumina/assembly/{sample}/contigs.gfa" )
    output:
        temp_assembly_dir=temp( directory( "intermediates/illumina/assembly/{sample}/" ) ),
        assembly="intermediates/illumina/assembly/{sample}_contigs.fasta",
        assembly_graph="intermediates/illumina/assembly/{sample}_contigs.gfa"
    threads: 8
    shell:
        """
        shovill \
            --assembler {params.assembler} \
            {params.assembler_arguments} \
            --outdir {output.temp_assembly_dir} \
            --force \
            --R1 {input.read1} \
            --R2 {input.read2} \
            --minlen {params.minimum_contig_length} \
            {params.trim} \
            {params.disable_correction} \
            {params.disable_read_correction} \
            {params.disable_read_stitching} &&\
        mv {params.temp_contigs} {output.assembly} &&\
        mv {params.temp_graph} {output.assembly_graph}
        """


rule contig_assignment:
    message: "Identify features of interest in contigs from {wildcards.sample}"
    input:
        contigs=rules.denovo_assembly.output.assembly
    params:
        prefix="intermediates/illumina/assignment/{sample}",
        temp_annotations="intermediates/illumina/assignment/{sample}.gff",
        prokka_arguments=config["contig_assignment"]["prokka_arguments"]
    output:
        annotated_contigs="results/assembly/{sample}.annotated.gff",
        prokka_output=expand(
            "intermediates/illumina/assignment/{{sample}}.{extension}",extension=["gbk", "fna", "faa",
                                                                                  "ffn", "sqn", "fsa",
                                                                                  "tbl", "err", "log",
                                                                                  "txt", "tsv"]
        )
    threads: min( workflow.cores,8 )
    shell:
        """
        prokka \
            --prefix {params.prefix} \
            --cpus {threads} \
            {input.contigs} &&\
        mv {params.temp_annotations} {output.annotated_contigs}
        """


include: "quality-control.smk"
