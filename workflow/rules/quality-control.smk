rule fastqc:
    input:
        read1s=[file["read1"] for file in SAMPLES.values()],
        read2s=[file["read2"] for file in SAMPLES.values()]
    output:
        directory=directory( "results/reports/fastqc/{sample}/" ),
    threads: 8
    shell:
        """
        fastqc \
            --outdir {output.directory} \
            --quiet \
            {input.read1s} {input.read2s}
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

#rule bamqc:
#    input:
#        alignment = rules.alignment_bwa.output.alignment
#    output:


rule combine_reports:
    input:
        expand( "results/reports/fastqc/{sample}/",sample=SAMPLES ),
        expand( "results/reports/samtools/{sample}.stats.txt",sample=SAMPLES ),
        expand( "results/reports/samtools/{sample}.idxstats.txt",sample=SAMPLES )
    output:
        report="results/reports/qc_report.html",
        report_directory=directory( "results/reports/qc_data/" )
    shell:
        """
        multiqc \
            --filename {output.report} \
            --outdir {output.report_directory} \
            results/reports/
        """
