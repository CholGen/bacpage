rule combine_index_genes:
    message: "Combine reference genes for typing into a single indexed file."
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
    message: "Align reads from {wildcards.sample} to typing reference genes."
    input:
        reference=rules.combine_index_genes.output.core_genes,
        reference_index=rules.combine_index_genes.output.index,
        reads1=lambda wildcards: SAMPLES[wildcards.sample]["read1"],
        reads2=lambda wildcards: SAMPLES[wildcards.sample]["read2"]
    params:
        bwa_params=config["alignment_bwa"]["bwa_params"]
    output:
        alignment=temp( "intermediates/illumina/typing/{sample}.typing.bam" ),
        alignment_index=temp( "intermediates/illumina/typing/{sample}.typing.bam.bai" )
    threads: min( 8,workflow.cores )
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
    message: "Determine which typing genes have coverage in {wildcards.sample}."
    input:
        alignment=rules.align_to_genes.output.alignment
    params:
        script_location=os.path.join( workflow.basedir,"scripts/calculate_typing_stats.py" ),
        minimum_coverage=config["coverage_mask"]["required_depth"]
    output:
        stats="intermediates/illumina/typing_stats/{sample}.stats.csv"
    shell:
        """
        python {params.script_location} \
            --alignment {input.alignment} \
            --min-coverage {params.minimum_coverage} \
            --sample-name {wildcards.sample} \
            --output {output.stats}
        """


rule virulence_factor_profiling:
    message: "Combine the virulence factor stats for each sample into a single report."
    input:
        stats=expand( "intermediates/illumina/typing_stats/{sample}.stats.csv",sample=SAMPLES )
    output:
        report="results/reports/typing_information.csv",
        detailed_report="results/reports/typing_information_detailed.csv"
    run:
        import pandas as pd

        stats = [pd.read_csv( file ) for file in input.stats]
        stats = pd.concat( stats )
        stats.to_csv( output.detailed_report,index=False )
        stats = stats.pivot( index="sample",columns="gene",values="frac_covered" ).sort_index()
        stats.to_csv( output.report )


rule mlst_profiling:
    message: "Scan consensus sequences against PubMLST typing schemes to identify sequence types."
    input:
        sequences=expand( "intermediates/illumina/consensus/{sample}.consensus.fasta",sample=SAMPLES )
    params:
        scheme=config["mlst_profiling"]["scheme"],
        mlst_params=config["mlst_profiling"]["mlst_params"]
    output:
        types="results/reports/mlst_types.csv"
    shell:
        """
        mlst \
            --scheme {params.scheme} \
            {params.mlst_params} \
            {input.sequences} > {output.types}
        """

rule antibiotic_resistance_profiling:
    input:
        consensus_sequences = expand( "results/consensus/{sample}.consensus.fasta",sample=SAMPLES )
    params:
        db = config["antibiotic_resistance"]["database"]
    output:
        temp_report = temp( "results/reports/antibiotic_resistance.temp.tab" ),
        summary = "results/reports/antibiotic_resistance.tsv"
    threads: min( workflow.cores, 8)
    shell:
        """
        abricate {input.consensus_sequences} \
            --threads {threads} \
            --nopath \
            --db {params.db} > {output.temp_report} && \
        abricate --summary {output.temp_report} > {output.summary}
        """
