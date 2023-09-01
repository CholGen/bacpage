def calculate_complete_sequences( wildcards ):
    depths = pd.read_csv( checkpoints.combine_depth.get( **wildcards ).output.combined )
    depths = depths.loc[depths["frac_covered"] > config["tree_building"]["minimum_completeness"]]
    complete_sequences = depths["sample"].apply( lambda x: f"results/consensus_sequences/{x}.masked.fasta" ).to_list()
    return complete_sequences


rule concatenate_sequences:
    input:
        background=config["background_dataset"],
        complete_sequences=calculate_complete_sequences
    output:
        alignment=temp( "intermediates/illumina/phylogeny/complete_alignment.fasta" )
    shell:
        """
        cat {input.background} {input.complete_sequences} > {output.alignment}
        """

rule reduce_alignment:
    message: "Removes redundant columns from the alignment, in order to save space reduce computation cost."
    input:
        alignment=rules.concatenate_sequences.output.alignment
    output:
        reduced_alignment="intermediates/illumina/phylogeny/reduced_alignment.fasta"
    shell:
        """
        snp-sites -m -o {output.reduced_alignment} {input.alignment}
        """

# TODO: Add substitution model. GTR+G is probably fine.
rule generate_tree:
    input:
        alignment=rules.reduce_alignment.output.reduced_alignment
    params:
        outgroup=config["tree_building"]["outgroup"],
        iqtree_parameters=config["tree_building"]["iqtree_parameters"]
    output:
        tree="intermediates/illumina/phylogeny/complete_alignment.fasta.treefile"
    threads: min( 16,workflow.cores )
    shell:
        """
        iqtree \
            -s {input.alignment} \
            -o {params.outgroup:q} \
            {params.iqtree_parameters}
        """
