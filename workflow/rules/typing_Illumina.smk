import os

rule index_genes:
    input:
        sequence=lambda wildcards: GENES[wildcards.gene]
    output:
        index=os.path.join( config["reference_genes"],"{gene}.fasta.ann" )
    shell:
        """
        bwa index {input.sequence}
        """

rule align_to_genes:
    input:
        reference=lambda wildcards: GENES[wildcards.gene],
        reference_index=os.path.join( config["reference_genes"],"{gene}.fasta.ann" ),
        reads1=lambda wildcards: SAMPLES[wildcards.sample]["read1"],
        reads2=lambda wildcards: SAMPLES[wildcards.sample]["read2"]
    params:
        bwa_params=config["alignment_bwa"]["bwa_params"]
    output:
        alignment="intermediates/illumina/typing/{sample}.{gene}.bam"
    threads: 8
    shell:
        """
        bwa mem \
            {params.bwa_params} \
            -t {threads} \
            {input.reference} \
            {input.reads1} \
            {input.reads2} | \
        samtools view -Sb - | \
        samtools fixmate - - | \
        samtools sort - > {output.alignment}
        """

rule calculate_depth:
    input:
        alignment=rules.align_to_genes.output.alignment
    output:
        depth="intermediates/illumina/typing_depth/{sample}.{gene}.depth.txt"
    shell:
        """
        samtools depth \
            -a {input.alignment} > {output.depth}
        """

rule calculate_gene_alignment_states:
    input:
        alignment=rules.align_to_genes.output.alignment,
        depth=rules.calculate_depth.output.depth
    params:
        minimum_coverage=config["coverage_mask"]["required_depth"]
    output:
        stats="intermediates/illumina/typing_stats/{sample}.{gene}.stats.csv"
    shell:
        """
        reads_mapped=$(samtools view -c -F 4 {input.alignment}) && \
        reads=$(samtools view -c {input.alignment}) && \
        frac_mapped=$(echo "scale=7; $reads_mapped / $reads" | bc) && \
        mean_depth=$(awk '{{ total += $3 }} END {{ print total/NR }}' {input.depth}) && \
        pos_covered=$(awk -v MIN_COV={params.minimum_coverage} 'BEGIN {{count = 0}} $3 > MIN_COV {{count++}} END {{print count}}' {input.depth}) && \
		total_pos=$(wc -l < {input.depth}) && \
		frac_covered=$(echo "scale=7; $pos_covered / $total_pos" | bc) && \
		printf "SAMPLENAME\tgene\treads_mapped\treads\tfrac_mapped\tmean_depth\tpos_covered\ttotal_pos\tfrac_covered\n" > {output.stats} && \
	    printf "{wildcards.sample}\t{wildcards.gene}\t$reads_mapped\t$reads\t$frac_mapped\t$mean_depth\t$pos_covered\t$total_pos\t$frac_covered\n" >> {output.stats}
        """


rule combine_gene_alignment_stats:
    input:
        stats=expand( "intermediates/illumina/typing_stats/{{sample}}.{gene}.stats.csv",gene=GENES )
    output:
        report="results/reports/{sample}.typing.csv"
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input.stats} > {output.report}
        """
