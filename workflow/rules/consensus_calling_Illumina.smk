import os.path

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


# Making the output temporary can save space but prevents pipeline from being run piecewise.
rule alignment_bwa:
    message: "Mapping reads for {wildcards.sample} to {input.reference} using `bwa mem`."
    input:
        reads1=lambda wildcards: SAMPLES[wildcards.sample]["read1"],
        reads2=lambda wildcards: SAMPLES[wildcards.sample]["read2"],
        reference=config["reference"],
        reference_index=rules.index_reference.output.reference_index
    params:
        bwa_params=config["alignment_bwa"]["bwa_params"]
    output:
        alignment="intermediates/illumina/merged_aligned_bams/{sample}.sorted.bam"
    threads: min( 8,workflow.cores )
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
            -o {output.alignment} - 
        """

# TODO: bedtools merge can be used to make depth_mask smaller.
rule generate_low_coverage_mask:
    message: "Create bed file from bam file for {wildcards.sample} indicating sites covered by less than {params.minimum_depth} reads"
    input:
        alignment=rules.alignment_bwa.output.alignment
    output:
        depth="intermediates/illumina/depth/{sample}.depth",
        depth_mask="intermediates/illumina/depth/{sample}.depthmask.bed"
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
            - > {output.depth_mask}
        """


checkpoint combine_depth:
    message: "Calculate percentage of reference thats covered by each sample"
    input:
        depth=expand( "intermediates/illumina/depth/{sample}.depth",sample=SAMPLES )
    params:
        minimum_depth=config["coverage_mask"]["required_depth"]
    output:
        combined="intermediates/misc/coverage.csv"
    run:
        import pandas as pd
        import os

        df = list()
        for file in input.depth:
            temp = pd.read_csv( file,sep="\t",names=["genome", "position", "depth"] )
            temp["sample"] = os.path.split( file )[-1].replace( ".depth","" )
            df.append( temp )
        df = pd.concat( df )
        df = df.groupby( "sample" )["depth"].agg( ["count", "mean", lambda x: sum( x > params.minimum_depth )] )
        df.columns = ["total_pos", "mean_depth", "pos_covered"]
        df["frac_covered"] = df["pos_covered"] / df["total_pos"]
        df.to_csv( output.combined )


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
    threads: min( 8,workflow.cores )
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
        variants=rules.call_variants_from_alignment.output.variants,
    params:
        minimum_depth=config["filter_variants"]["minimum_depth"],
        minimum_strand_depth=config["filter_variants"]["minimum_strand_depth"],
        minimum_support=config["filter_variants"]["minimum_support"]
    output:
        filtered_variants="intermediates/illumina/variants/{sample}.bcftools.filt.vcf"
    shell:
        """
        bcftools filter \
            --no-version \
            -i "INFO/AD[1]>{params.minimum_depth} && (INFO/AD[1])/(INFO/AD[0]+INFO/AD[1])>{params.minimum_support} && INFO/ADF[1]>{params.minimum_strand_depth} && INFO/ADR[1]>{params.minimum_strand_depth}" \
            -o {output.filtered_variants} \
            {input.variants}
        """

rule align_and_normalize_variants:
    message: "For sample {wildcards.sample}, Left-align and normalize indels, and remove insertions."
    input:
        variants=rules.filter_variants.output.filtered_variants,
        reference=config["reference"],
        reference_index=config["reference"] + ".bwt"
    output:
        normalized_variants="intermediates/illumina/variants/{sample}.bcftools.filt.norm.vcf.gz",
        variant_index="intermediates/illumina/variants/{sample}.bcftools.filt.norm.vcf.gz.csi"
    shell:
        """
        bcftools norm \
            --no-version \
            -f {input.reference} \
            {input.variants} |\
        bcftools filter \
            --no-version \
            --exclude 'strlen(REF)<strlen(ALT)' \
            -Oz \
            -o {output.normalized_variants} &&\
        bcftools index {output.normalized_variants}
        """


rule call_consensus:
    message: "For sample {wildcards.sample}, apply variants to reference to create consensus sequences. Masks sites with less than desired coverage."
    input:
        variants=rules.align_and_normalize_variants.output.normalized_variants,
        variant_index=rules.align_and_normalize_variants.output.variant_index,
        depth_mask=rules.generate_low_coverage_mask.output.depth_mask,
        reference=config["reference"],
        reference_index=config["reference"] + ".bwt"
    params:
        bcftools_parameters=config["call_consensus"]["consensus_parameters"]
    output:
        consensus_sequence="results/consensus/{sample}.consensus.fasta"
    shell:
        """
	    bcftools consensus \
	        {params.bcftools_parameters} \
	        --fasta-ref {input.reference} \
	        --mask {input.depth_mask} \
	        {input.variants} |\
        union -filter |\
        sed "1s/.*/> {wildcards.sample}/" > {output.consensus_sequence}
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
