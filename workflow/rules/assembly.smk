# Commented out rules are not implemented yet.

rule basecalling:
    input:
        run_configuration = config["path-to-run-configuration"]
    params:
        run_directory = config["path-to-raw-data"],
        guppy_params = config["basecalling"]["guppy_params"]
    output:
        output_directory = directory( "input/fastqs" ),
        summary = "input/fastqs/sequencing_summary.txt"
    log: "results/logs/basecalling.txt"
    shell:
        """
        guppy_basecaller \ 
            {params.guppy_params} \
            -c {input.run_configuration} \ 
            -i {params.run_directory} \
            -s {output.output_directory} 2>&1 {log}
        """


rule visualize_run_metrics:
    input:
        basecalls = rules.basecalling.output.output_directory,
        sequencing_summary = rules.basecalling.output.summary
    output:
        summary_directory = "results/reports/basecalling",
        summary_html = "results/reports/basecalling/<x>.html"
    shell:
        """
        NanoPlot \
            --summary {input.sequencing_summary} \
            -o {output.summary_directory}
        """


# rule demultiplex:


rule filter_reads:
    message: "Filter reads from {wildcards.sample} if filesizes are too large. Keep {params.keep_percent}% of reads, up to {params.target_bases:,d} bases."
    input:
        fastq = "input/fastqs/{sample}_1.fastq.gz"
    params:
        keep_percent = config["filter_reads"]["keep_percent"],
        target_bases = config["filter_reads"]["target_bases"],
        extra_params = config["filter_reads"]["filtlong_params"]
    output:
        fastq = "intermediates/filter_reads/{sample}.filtered.fastq.gz"
    shell:
        """
        filtlong \
            --keep_percent {params.keep_percent} \
            --target_bases {params.target_bases} \
            {params.extra_params} \
            {input.fastq} | \
        gzip > {output.fastq}
        """


# TODO: this seems like it should be done, but I don't have access to the run directory at the moment.
rule index_reads:
    message: "To aid alignment, index reads for sample {wildcards.sample}."
    input:
        raw_reads = config["path-to-raw-data"],
        read_summary = rules.basecalling.output.summary,
        reads = rules.filter_reads.output.fastq
    output:
        index = "intermediates/filter_reads/{sample}.filtered.fastq.gz.index"
    shell:
        """
        nanopolish index \
            -d {input.raw_reads} \
            -s {input.read_summary} \
            {input.reads}
        """


rule alignment:
    message: "Align reads in {wildcards.sample} to reference: {input.reference}"
    input:
        reads = rules.filter_reads.output.fastq,
        reference = config["reference"]
    output:
        alignment = "intermediates/alignment/{sample}.sorted.bam",
        index = "intermediates/alignment/{sample}.sorted.bam.bai"
    threads: 16
    shell:
        """
        minimap2 -a -x map-ont -t {threads} {input.reference} {input.reads} | \
        samtools sort \
            --threads {threads} \
            -o {output.alignment} && \
        samtools index {output.alignment}
        """
