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


rule convert_gff_to_bed:
    input:
        gff=config["recombinant_mask"]
    output:
        bed="intermediates/misc/recombinant_mask.bed"
    run:
        import pandas as pd

        gff = pd.read_csv( input.gff,sep="\t",header=None )
        bed = gff[[0, 3, 4]].copy()
        bed[3] -= 1
        bed.to_csv( output.bed,sep="\t",header=False,index=False )


rule convert_alignment_to_vcf:
    message: "Converts multiple sequence alignment to sparse alignment, in VCF format."
    input:
        alignment=rules.concatenate_sequences.output.alignment
    params:
        reference=config["tree_building"]["outgroup"]
    output:
        vcf=temp( "intermediates/illumina/phylogeny/alignment.vcf" )
    shell:
        """
        faToVcf \
            -includeRef -ambiguousToN \
            -ref={params.reference:q} \
            {input.alignment} {output.vcf}
        """


rule mask_vcf:
    message: "Masking estimated recombinant sites from all sequences in alignment."
    input:
        vcf=rules.convert_alignment_to_vcf.output.vcf,
        mask=rules.convert_gff_to_bed.output.bed
    output:
        masked_vcf="intermediates/illumina/phylogeny/alignment.masked.vcf"
    shell:
        """
        augur mask \
            --sequences {input.vcf} \
            --mask {input.mask} \
            --output {output.masked_vcf}
        """


rule concatenate_reference:
    input:
        reference=config["reference"]
    output:
        reference="intermediates/misc/concatenated_reference.fasta"
    shell:
        """
        union -filter {input.reference} > {output.reference}
        """


rule generate_tree:
    input:
        alignment=rules.mask_vcf.output.masked_vcf,
        reference=rules.concatenate_reference.output.reference
    params:
        iqtree_parameters=config["tree_building"]["iqtree_parameters"]
    output:
        tree="intermediates/misc/initial_phylogeny.tree"
    threads: min( 16,workflow.cores )
    shell:
        """
        augur tree \
            --substitution-model GTR \
            --nthreads {threads} \
            --alignment {input.alignment} \
            --vcf-reference {input.reference} \
            --tree-builder-args={params.iqtree_parameters:q} \
            --override-default-args \
            --output {output.tree}
        """


rule generate_rooted_tree:
    input:
        tree=rules.generate_tree.output.tree
    params:
        outgroup=config["tree_building"]["outgroup"]
    output:
        rooted_tree="results/phylogeny/phylogeny.tree"
    run:
        from Bio import Phylo

        tree = Phylo.read( input.tree,"newick" )
        tree.root_with_outgroup( params.outgroup )
        Phylo.write( tree,output.rooted_tree,"newick" )
