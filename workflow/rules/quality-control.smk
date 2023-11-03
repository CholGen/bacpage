rule rename_fastq:
    message: "Generate temporary copy of FASTQ file to maintain names for QC."
    input:
        reads1=lambda wildcards: config["SAMPLES"][wildcards.sample]["read1"]
    output:
        renamed_reads=temp( "intermediates/tmp/{sample}.fastq.gz" )
    shell:
        """
        cp {input.reads1} {output.renamed_reads}
        """


rule fastqc:
    # We're assuming that R1 is representative of R2. This should generally work and I can't think of a reason where
    # problems would only pop up in one rather than the other.
    message: "Calculate quality control metrics for raw sequencing reads of {wildcards.sample}."
    input:
        reads1=rules.rename_fastq.output.renamed_reads
    output:
        directory=directory( "results/reports/fastqc/{sample}/" ),
    threads: min( 2,workflow.cores )
    shell:
        """
        mkdir {output.directory} && \
        fastqc \
            --outdir {output.directory} \
            --threads {threads} \
            --quiet \
            {input.reads1}  
        """


rule alignment_stats:
    input:
        alignment="intermediates/illumina/merged_aligned_bams/{sample}.sorted.bam"
    output:
        alignment_stats="results/reports/samtools/{sample}.stats.txt",
        alignment_idxstats="results/reports/samtools/{sample}.idxstats.txt"
    shell:
        """
        samtools index {input.alignment} && \
        samtools idxstats {input.alignment} > {output.alignment_idxstats} && \
        samtools stats {input.alignment} > {output.alignment_stats} 
        """


rule bamqc:
    input:
        alignment="intermediates/illumina/merged_aligned_bams/{sample}.sorted.bam"
    output:
        reheaded_alignment="intermediates/illumina/merged_aligned_bams/{sample}.headed.bam",
        report_directory=directory( "results/reports/bamqc/{sample}/" )
    threads: min( 8,workflow.cores )
    shell:
        """
        samtools view -H {input.alignment} |\
        sed 's,^@RG.*,@RG\\tID:None\\tSM:None\\tLB:None\\tPL:Illumina,g' |\
        samtools reheader - {input.alignment} > {output.reheaded_alignment} && \
        qualimap bamqc \
            -bam {output.reheaded_alignment} \
            -nt {threads} \
            -outdir {output.report_directory}
        """


rule assembly_stats:
    input:
        assembly="intermediates/illumina/assembly/{sample}.assembly.fasta"
    params:
        quast_arguments="--fast --space-efficient"
    output:
        report_directory=directory( "results/reports/quast/{sample}/" )
    threads: min( 8,workflow.cores )
    shell:
        """
        quast \
            {params.quast_arguments} \
            -o {output.report_directory} \
            {input.assembly}
        """


rule coverage_plot:
    input:
        depth="intermediates/illumina/depth/{sample}.depth"
    params:
        script_location=workflow.source_path( "../scripts/plot_coverage.py" ),
        bin_size=config["plot_coverage"]["bin_size"],
        minimum_depth=config["coverage_mask"]["required_depth"]
    output:
        coverage_plot="results/reports/depth/{sample}.depth.pdf"
    shell:
        """
        python {params.script_location} \
            --input {input.depth} \
            --bin-size {params.bin_size} \
            --min-depth {params.minimum_depth} \
            --output {output.coverage_plot}
        """


def get_qc_inputs( wildcards ):
    inputs = list()
    if config["DENOVO"]:
        inputs.extend( expand( "results/reports/fastqc/{sample}/",sample=config["SAMPLES"] ) )
        inputs.extend( expand( "results/reports/quast/{sample}/",sample=config["SAMPLES"] ) )
    else:
        inputs.extend( expand( "results/reports/fastqc/{sample}/",sample=config["SAMPLES"] ) )
        inputs.extend( expand( "results/reports/samtools/{sample}.stats.txt",sample=config["SAMPLES"] ) )
        inputs.extend( expand( "results/reports/samtools/{sample}.idxstats.txt",sample=config["SAMPLES"] ) )
        inputs.extend( expand( "results/reports/bamqc/{sample}/",sample=config["SAMPLES"] ) )
        inputs.extend( expand( "results/reports/depth/{sample}.depth.pdf",sample=config["SAMPLES"] ) )
    return inputs


rule generate_complete_report:
    input:
        get_qc_inputs
    params:
        multiqc_config=workflow.source_path( "../../resources/multiqc_config.yaml" )
    output:
        report="results/reports/qc_report.html",
        report_directory=directory( "results/reports/qc_report_data/" )
    shell:
        """
        multiqc \
            --filename {output.report} \
            --config {params.multiqc_config} \
            results/reports/
        """
