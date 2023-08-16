rule combine_index_genes:
    input:
        sequence=GENES.values()
    output:
        core_genes="intermediates/illumina/typing/reference_genes.fasta",
        index=multiext( "intermediates/illumina/typing/reference_genes.fasta",".bwt",".pac",".ann",".sa",".amb" )
    run:
        from Bio import SeqIO

        all_seqs = list()
        for name, seq in GENES.items():
            record = SeqIO.read( seq,"fasta" )
            record.id = name
            record.name = ""
            record.description = ""
            all_seqs.append( record )
        SeqIO.write( all_seqs,output.core_genes,"fasta" )

        shell( "bwa index {output.core_genes}" )

rule align_to_genes:
    input:
        reference=rules.combine_index_genes.output.core_genes,
        reference_index=rules.combine_index_genes.output.index,
        reads1=lambda wildcards: SAMPLES[wildcards.sample]["read1"],
        reads2=lambda wildcards: SAMPLES[wildcards.sample]["read2"]
    params:
        bwa_params=config["alignment_bwa"]["bwa_params"]
    output:
        alignment="intermediates/illumina/typing/{sample}.typing.bam",
        alignment_index="intermediates/illumina/typing/{sample}.typing.bam.bai"
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
        samtools sort - > {output.alignment} && \
        samtools index {output.alignment}
        """


rule calculate_gene_alignment_states:
    input:
        alignment=rules.align_to_genes.output.alignment
    params:
        minimum_coverage=config["coverage_mask"]["required_depth"]
    output:
        stats="intermediates/illumina/typing_stats/{sample}.stats.csv"
    shell:
        """
        python workflow/scripts/calculate_typing_stats.py \
            --alignment {input.alignment} \
            --min-coverage {params.minimum_coverage} \
            --sample-name {wildcards.sample} \
            --output {output.stats}
        """

rule combine_gene_alignment_stats:
    input:
        stats=expand( "intermediates/illumina/typing_stats/{sample}.stats.csv",sample=SAMPLES )
    output:
        report="results/reports/typing_information.csv"
    shell:
        """
        awk '(NR == 1) || (FNR > 1)' {input.stats} > {output.report}
        """
