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
    params:
        minimmap_params = config["alignment"]["minimap_params"]
    output:
        alignment = "intermediates/alignment/{sample}.sorted.bam",
        index = "intermediates/alignment/{sample}.sorted.bam.bai"
    threads: 16
    shell:
        """
        minimap2 -a -x map-ont {params.minimmap_params} -t {threads} {input.reference} {input.reads} | \
        samtools sort \
            --threads {threads} \
            -o {output.alignment} && \
        samtools index {output.alignment}
        """


rule calculate_depth:
    message: "Calculate the read depth at each position of the alignment for {wildcards.sample}"
    input:
        alignment = rules.alignment.output.alignment
    params:
        maximum_depth = config["calculate_depth"]["maximum_depth"]
    output:
        depth = "intermediates/depth/{sample}.depth.txt"
    shell:
        """
        samtools depth \
            -a -m {params.maximum_depth} \
            {input.alignment} > {output.depth}
        """


checkpoint calculate_windows:
    message: "Break up reference into window sizes to prevent loading all fast5 data into memory each time."
    input:
        reference = config["reference"]
    output:
        windows = "intermediates/ref_windows.txt"
    shell:
        """
        nanopolish_makerange.py {input.reference} > {output.windows}
        """


def get_window( wildcards ):
    windows_file = checkpoints.calculate_windows.get( **wildcards ).output[0]
    with open( windows_file, "r" ) as windows:
        for idx, line in enumerate( windows ):
            if idx == i:
                return line.strip()

# Cannot run this rule until I have sequencing summary.
rule nanopolish:
    message: "Calculate SNPs using a signal-level HMM for window {params.window} of {wildcards.sample}"
    input:
        windows = rules.calculate_windows.output.windows,
        reads = rules.filter_reads.output.fastq,
        index = rules.index_reads.output.index,
        alignment = rules.alignment.output.alignment,
        reference = config["reference"]
    params:
        window = get_window
    output:
        vcf = temporary( "intermediates/nanopolish/{sample}.{count}.vcf" )
    threads: 1
    shell:
        """
        nanopolish variants \
            -o {output.vcf} \
            -p 1 \
            -w {params.window} \
            -q dam \
            -r {input.reads} \
            -b {input.alignment} \
            -g {input.reference} \
            -t {threads}
        """

def get_windows( wildcards ):
    windows_file = checkpoints.calculate_windows.get( **wildcards ).output[0]
    windows = []
    with open( windows_file, "r" ) as wfile:
        for idx, _ in enumerate( wfile ):
            windows.append( idx )
    return expand( "intermediates/nanopolish/{{sample}}.{count}.vcf", count = windows )


rule combine_vcfs:
    message: "Concatenate regional VCF files for sample: {wildcards.sample}"
    input:
        reference = config["reference"],
        vcfs = get_windows
    output:
        concated_vcf = "intermediates/variants/{sample}.vcf"
    shell:
        """
        python concatenate_vcfs.py \
            --reference {input.reference} \
            --vcfs {input.vcfs} \
            --output {output.concated_vcf}
        """

rule filter_vcf:
    message: "Filter variants from VCF with depth less than {params.minimum_reads}, and support less than {params.minimum_support:.0%}"
    input:
        vcf = rules.combine_vcfs.output.concated_vcf
    params:
        minimum_reads = config["filter_vcf"]["minimum_reads"],
        minimum_support = config["filter_vcf"]["minimum_supports"]
    output:
        filtered_vcf = temporary( "intermediates/variants/{sample}.filtered.vcf.gz" )
    shell:
        """
        bcftools filter \
            --no-version \
            -i "INFO/SupportFraction>{params.minimum_support}" \ 
            {input.vcf} | \
		bcftools filter \
		    --no-version \
		    -i "INFO/TotalReads>{params.minimum_reads}" \
		    -O z \
		    -o {output.filtered_vcf}
        """

rule rename_index_vcf:
    message: "Rename and index VCF file for later manipulation"
    input:
        vcf = rules.filter_vcf.output.filtered_vcf
    output:
        header = temporary( "intermediates/variants/{sample}.header.txt" ),
        renamed_vcf = "intermediates/variants/{sample}.renamed.vcf.gz"
    threads: 1
    shell:
        """
        echo "sample {wildcards.sample}" > {output.header} &&
        bcftools reheader \
            -s {output.header} \
            -o {output.renamed_vcf} \
            {input.vcf} &&
        bcftools index \
            --threads {threads} \
            {output.renamed_vcf}
        """
