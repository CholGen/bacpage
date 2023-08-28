rule fastqc:
    # We're assuming that R1 is representative of R2. This should generally work and I can't think of a reason where
    # problems would only pop up in one rather than the other.
    message: "Calculate quality control metrics for raw sequencing reads of {wildcards.sample}."
    input:
        reads1=lambda wildcards: SAMPLES[wildcards.sample]["read1"],
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
        alignment=rules.alignment_bwa.output.alignment
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
        alignment=rules.alignment_bwa.output.alignment
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


rule generate_complete_report:
    input:
        expand( "results/reports/fastqc/{sample}/",sample=SAMPLES ),
        expand( "results/reports/samtools/{sample}.stats.txt",sample=SAMPLES ),
        expand( "results/reports/samtools/{sample}.idxstats.txt",sample=SAMPLES ),
        expand( "results/reports/bamqc/{sample}/",sample=SAMPLES )
    params:
        multiqc_config=os.path.join( workflow.basedir,"../resources/multiqc_config.yaml" )
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
