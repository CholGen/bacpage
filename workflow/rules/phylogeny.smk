rule all:
    input:
        tree="results/phylogeny/phylogeny.tree"


rule mask_consensus:
    message: "Mask wholly recombinant regions of the consensus sequence of {wildcards.sample} using the mask provided in {input.recombinant_mask}"
    input:
        consensus=lambda wildcards: config["SAMPLES"][wildcards.sample],
        recombinant_mask=config["recombinant_mask"]
    output:
        first_tmp_consensus=temp( "intermediates/illumina/consensus/{sample}.tmp1.fasta" ),
        second_tmp_consensus=temp( "intermediates/illumina/consensus/{sample}.tmp2.fasta" ),
        masked_consensus="results/consensus/{sample}.masked.fasta"
    shell:
        """
        HEADER=$(cut -f1 {input.recombinant_mask} | uniq | head -n1) && \
        sed "1s/.*/>${{HEADER}}/" {input.consensus} > {output.first_tmp_consensus} && \
        bedtools maskfasta \
            -fi {output.first_tmp_consensus} \
            -bed {input.recombinant_mask} \
            -fo {output.second_tmp_consensus} && \
        sed "1s/.*/>{wildcards.sample}/" {output.second_tmp_consensus} > {output.masked_consensus}
    """


def calculate_complete_sequences( wildcards ):
    complete_sequences = expand( "results/consensus/{sample}.masked.fasta",sample=config["SAMPLES"] )
    if config["BACKGROUND"] != "":
        complete_sequences.append( config["BACKGROUND"] )
    return complete_sequences


rule concatenate_sequences:
    input:
        sequences=calculate_complete_sequences
    output:
        alignment=temp( "intermediates/illumina/phylogeny/complete_alignment.fasta" )
    shell:
        """
        cat {input.sequences} > {output.alignment}
        """


rule generate_sparse_alignment:
    input:
        alignment=rules.concatenate_sequences.output.alignment
    output:
        sparse_alignment="results/phylogeny/sparse_alignment.fasta"
    shell:
        """
        snp-sites -o {output.sparse_alignment} {input.alignment}
        """


def calculate_outgroup( wildcards ):
    outgroup = config["generate_tree"]["outgroup"]
    return f"-o {outgroup:q}" if outgroup != "" else ""


rule generate_tree:
    input:
        alignment=rules.generate_sparse_alignment.output.sparse_alignment,
    output:
        tree=temp(
            expand(
                "results/phylogeny/sparse_alignment.fasta" + '.{extension}',extension=["iqtree", "treefile", "mldist",
                                                                                       "splits.nex", "contree", "log"]
            )
        )
    threads: min( 16,workflow.cores )
    shell:
        """
        iqtree \
            -nt {threads} \
            -m GTR \
            -B 1000 \
            -s {input.alignment}
        """


rule move_tree_and_rename:
    input:
        iqtree_output="results/phylogeny/sparse_alignment.fasta.treefile"
    output:
        final_tree="results/phylogeny/phylogeny.tree"
    shell:
        """
        mv {input.iqtree_output} {output.final_tree}
        """

    #rule convert_gff_to_bed:
    #    input:
    #        gff=config["recombinant_mask"]
    #    output:
    #        bed="intermediates/misc/recombinant_mask.bed"
    #    run:
    #        import pandas as pd
    #
    #        gff = pd.read_csv( input.gff,sep="\t",header=None )
    #        bed = gff[[0, 3, 4]].copy()
    #        bed[3] -= 1
    #        bed.to_csv( output.bed,sep="\t",header=False,index=False )


    #rule convert_alignment_to_vcf:
    #    message: "Converts multiple sequence alignment to sparse alignment, in VCF format."
    #    input:
    #        alignment=rules.concatenate_sequences.output.alignment
    #    params:
    #        reference=config["tree_building"]["outgroup"],
    #        script_location = os.path.join( workflow.basedir,"scripts/faToVcf" )
    #    output:
    #        vcf=temp( "intermediates/illumina/phylogeny/alignment.vcf" )
    #    shell:
    #        """
    #        {params.script_location} \
    #            -includeRef -ambiguousToN \
    #            -ref={params.reference:q} \
    #            {input.alignment} {output.vcf}
    #        """

    #rule convert_alignment_to_vcf:
    #    message: "Converts multiple sequence alignment to sparse alignment, in VCF format."
    #    input:
    #        alignment=rules.concatenate_sequences.output.alignment
    #    params:
    #        reference=config["tree_building"]["outgroup"]
    #    output:
    #        vcf=temp( "intermediates/illumina/phylogeny/alignment.vcf" )
    #    shell:
    #        """
    #        faToVcf \
    #            -includeRef -ambiguousToN \
    #            -ref={params.reference:q} \
    #            {input.alignment} {output.vcf}
    #        """


    #rule mask_vcf:
    #    message: "Masking estimated recombinant sites from all sequences in alignment."
    #    input:
    #        vcf=rules.convert_alignment_to_vcf.output.vcf,
    #        mask=rules.convert_gff_to_bed.output.bed
    #    output:
    #        masked_vcf="intermediates/illumina/phylogeny/alignment.masked.vcf"
    #    shell:
    #        """
    #        augur mask \
    #            --sequences {input.vcf} \
    #            --mask {input.mask} \
    #            --output {output.masked_vcf}
    #        """


    #rule concatenate_reference:
    #    input:
    #        reference=config["reference"]
    #    output:
    #        reference="intermediates/misc/concatenated_reference.fasta"
    #    shell:
    #        """
    #        union -filter {input.reference} > {output.reference}
    #        """


    # Add a conditional, if root is specified, root the tree, otherwise, just copy to the results.
    # TODO: Rename tree to name of project directory.
    #rule generate_rooted_tree:
    #    input:
    #        tree=rules.generate_tree.output.tree
    #    params:
    #        outgroup=config["tree_building"]["outgroup"]
    #    output:
    #        rooted_tree="results/phylogeny/phylogeny.tree"
    #    run:
    #        from Bio import Phylo
    #
    #        tree = Phylo.read( input.tree,"newick" )
    #        tree.root_with_outgroup( params.outgroup )
    #        Phylo.write( tree,output.rooted_tree,"newick" )
