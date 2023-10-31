def estimate_output( wildcards ):
    output = list()
    for sample in config["SAMPLES"]:
        output.append( f"results/assembly/{sample}.annotated.gff" )
    if config["QC"]:
        output.append( "results/reports/qc_report.html" )
    return output


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
        window_size=config["trimming"]["window_size"],
        window_minimum_quality=config["trimming"]["minimum_quality"],
        minimum_length=config["trimming"]["minimum_length"],
        baseout="intermediates/trimmed/{sample}.fastq.gz"
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
        firstpass_reference="adapters",
        firstpass_args="ktrim=r k=23 mink=11 hdist=1 tpe tbo",
        secondpass_reference="phix",
        secondpass_args="k=31 hdist=1"
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
        minimum_contig_length=200,
        trim=False,
        no_read_correction=False,
        no_stitch=False,
        no_corr=False,
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
            --outdir {output.temp_assembly_dir} \
            --force \
            --R1 {input.read1} \
            --R2 {input.read2} \
            --minlen {params.minimum_contig_length} &&\
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
        prokka_arguments=config["assignment"]["prokka_arguments"]
    output:
        annotated_contigs="results/assembly/{sample}.annotated.gff",
        prokka_output=expand( "intermediates/illumina/assignment/{{sample}}.{extension}",extension=["gbk", "fna", "faa",
                                                                                                    "ffn", "sqn", "fsa",
                                                                                                    "tbl", "err", "log",
                                                                                                    "txt", "tsv"] )
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
