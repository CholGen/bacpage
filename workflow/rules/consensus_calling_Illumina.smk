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


rule generate_low_coverage_mask:
    message: "Create bed file from bam file for {wildcards.sample} indicating sites covered by less than {params.minimum_depth} reads"
    input:
        alignment=rules.alignment_bwa.output.alignment
    output:
        depth="intermediates/illumina/depth/{sample}.depth",
        bed_file="intermediates/illumina/depth/{sample}.depthmask.bed"
    params:
        minimum_depth=config["coverage_mask"]["required_depth"]
    shell:
        """
        samtools depth \
            -a {input.alignment} |\
        tee {output.depth} |\
        awk \
            -v depth="{params.minimum_depth}" \
            '$3 < depth {{printf "%s\\t%d\\t%d\\n", $1, $2 - 1, $2}}' \
            - > {output.bed_file}
        """


rule call_variants_from_alignment:
    message: "Call variants from alignment for {wildcards.sample} using bcftools."
    input:
        alignment=rules.alignment_bwa.output.alignment,
        reference=config["reference"],
        reference_index=config["reference"] + ".bwt"
    params:
        maximum_depth=config["call_variants"]["maximum_depth"],
        minimum_mapping_quality=config["call_variants"]["minimum_mapping_quality"],
        minimum_base_quality=config["call_variants"]["minimum_base_quality"],
        mpileup_parameters=config["call_variants"]["mpileup_parameters"],
        call_parameters=config["call_variants"]["call_parameters"]
    output:
        variants="intermediates/illumina/variants/{sample}.bcftools.vcf"
    threads: 8
    shell:
        """
        bcftools mpileup \
            --threads {threads} \
            -d {params.maximum_depth} \
            -q {params.minimum_mapping_quality} \
            -Q {params.minimum_base_quality} \
            {params.mpileup_parameters} \
            -f {input.reference} \
            {input.alignment} |\
        bcftools call \
            --threads {threads} \
            {params.call_parameters} \
            -o {output.variants}
        """


rule filter_variants:
    message:
        """Remove variants for sample {wildcards.sample} that:
            - have depth less than {params.minimum_depth}
            - have individual strand depth less than {params.minimum_strand_depth}
            - are present in less than {params.minimum_support:.0%} of reads
        """
    input:
        variants=rules.call_variants_from_alignment,
    params:
        minimum_depth=config["filter_variants"]["minimum_depth"],
        minimum_strand_depth=config["filter_variants"]["minimum_strand_depth"],
        minimum_support=config["filter_variants"]["minimum_support"]
    output:
        filtered_variants="intermediates/illumina/variants/{sample}.filt.vcf"
    shell:
        """
            bcftools filter \
                --no-version \
                -i "INFO/AD[1]>{params.minimum_depth} && (INFO/AD[1])/(INFO/AD[0]+INFO/AD[1])>{params.minimum_support} && INFO/ADF[1]>{params.minimum_strand_depth} && INFO/ADR[1]>{params.minimum_strand_depth}" \
                -o {output.filtered_variants}
                {input.variants}
        """
    #rule trim_alignment:
    #    message: "Trim reads in alignment which are not a high quality."
    #    input:
    #        alignment=rules.alignment_bwa.output.alignment
    #    params:
    #        minimum_length=config["trimming"]["minimum_length"],
    #        minimum_quality=config["trimming"]["minimum_quality"],
    #        window_length=config["trimming"]["window_length"]
    #    output:
    #        trimmed_alignment="intermediates/illumina/trimmed_bams/{sample}.trimmed.sorted.bam"
    #    log: "results/logs/{sample}.trimming.log"
    #    shell:
    #        """
    #        ivar trim \
    #            -i {input.alignment} \
    #            -m {params.minimum_length} \
    #            -q {params.minimum_quality} \
    #            -s {params.window_length} 2> {log} |\
    #        samtools sort \
    #            -o {output.trimmed_alignment} -
    #        """
