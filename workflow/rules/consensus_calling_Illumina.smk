rule index_reference:
    input:
        reference=config["reference"]
    output:
        reference_index=expand(
            config["reference"] + '.{extension}',extension=["bwt", "pac", "ann", "amb", "sa"]
        )
    shell:
        """
        bwa index {input.reference}
        """

rule alignment_bwa:
    message: "Mapping reads for {wildcards.sample} to {input.reference} using `bwa mem`."
    input:
        reads1=lambda wildcards: SAMPLES[wildcards.sample]["read1"],
        reads2=lambda wildcards: SAMPLES[wildcards.sample]["read2"],
        reference=config["reference"],
        reference_index=config["reference"] + ".bwt"
    params:
        bwa_params=config["alignment_bwa"]["bwa_params"]
    output:
        alignment="intermediates/illumina/merged_aligned_bams/{sample}.sorted.bam",
        alignment_stats="intermediates/illumina/merged_aligned_bams/{sample}.stats"
    threads: 8
    shell:
        """
        bwa mem \
            {params.bwa_params} \
            -t {threads} \
            {input.reference} \
            {input.reads1} {input.reads2} |\
        samtools view -Sb - |\
        samtools sort - |\
        samtools addreplacerg \
            -r "ID:{wildcards.sample}" \
            -o {output.alignment} - && \
        samtools stats {output.alignment} > {output.alignment_stats}
        """

rule trim_alignment:
    message: "Trim reads in alignment which are not a high quality."
    input:
        alignment=rules.alignment_bwa.output.alignment
    params:
        minimum_length=config["trimming"]["minimum_length"],
        minimum_quality=config["trimming"]["minimum_quality"],
        window_length=config["trimming"]["window_length"]
    output:
        trimmed_alignment="intermediates/illumina/trimmed_bams/{sample}.trimmed.sorted.bam"
    log: "results/logs/{sample}.trimming.log"
    shell:
        """
        ivar trim \
            -i {input.alignment} \
            -m {params.minimum_length} \
            -q {params.minimum_quality} \
            -s {params.window_length} 2> {log} |\
        samtools sort \
            -o {output.trimmed_alignment} - 
        """
