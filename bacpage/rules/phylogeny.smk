def determine_output( wildcards ):
    if config["TERRA"]:
        return determine_alignment_from_vcf_input( wildcards )
    return "results/phylogeny/phylogeny.tree"


rule all:
    input:
        determine_output


def calculate_complete_sequences( wildcards ):
    complete_sequences = config["SAMPLES"].values()
    if config["BACKGROUND"] != "":
        if not config["BACKGROUND"].endswith( (".vcf", ".vcf.gz", ".bcf", ".bcf.gz") ):
            complete_sequences.append( config["BACKGROUND"] )
    return complete_sequences


rule concatenate_sequences:
    input:
        sequences=calculate_complete_sequences
    output:
        alignment=temp( "intermediates/illumina/alignment/complete_alignment.fasta" )
    shell:
        """
        cat {input.sequences} > {output.alignment}
        """


# from https://stackoverflow.com/a/69473926. Need to test at some point.
rule concatenate_reference:
    input:
        reference=config["reference"]
    output:
        concatenated_reference=temp( "intermediates/illumina/reference.fasta" )
    shell:
        """
        sed '1h;/>/d;H;$!d;x;s/\n//2g' {input.reference} > {output.concatenated_reference}
        """


rule convert_to_vcf:
    input:
        alignment=rules.concatenate_sequences.output.alignment,
        reference=rules.concatenate_reference.output.concatenated_reference
    output:
        chromosome_name=temp( "intermediates/illumina/alignment/chromosome_name.txt" ),
        temp_alignment=temp( "intermediates/illumina/alignment/reference_alignment.fasta" ),
        vcf="intermediates/illumina/alignment/complete_alignment.bcf.gz",
        vcf_index="intermediates/illumina/alignment/complete_alignment.bcf.gz.csi"
    shell:
        """
        REFERENCE=$(head -n1 {input.reference} | cut -f2 -d \> | cut -f2 -d" ") &&\
        echo "1 ${{REFERENCE}}" > {output.chromosome_name} &&\
        cat {input.reference} {input.alignment} > {output.temp_alignment} &&\
        snp-sites -v {output.temp_alignment} | bcftools annotate --rename-chrs {output.chromosome_name} -O b -o {output.vcf} &&\
        bcftools index {output.vcf}
        """


# TODO: Would love to test that this works as expected.
# TODO: test that phylogeny.py throws an error if given a vcf file as input and the #CHROM field doesn't match the reference being used.
rule combine_sequences_and_background_vcf:
    input:
        user_sequences=rules.convert_to_vcf.output.vcf,
        background_sequences=config["BACKGROUND"]
    output:
        combined_vcf="intermediates/illumina/alignment/combined_alignment.bcf.gz"
    shell:
        """
        bcftools merge \
            --output {output.combined_vcf} \
            --output-type b \
            {input.background_sequences} {input.user_sequences} 
        """


def determine_mask_vcf_inputs( wildcards ):
    if (config["BACKGROUND"] == "") or config["BACKGROUND"].endswith( (".fa", ".fasta") ):
        return rules.convert_to_vcf.output.vcf
    return rules.combine_sequences_and_background_vcf.output.combined_vcf


rule mask_vcf:
    input:
        alignment=determine_mask_vcf_inputs,
        mask=config["MASK_LOCATION"]
    output:
        masked_alignment="intermediates/illumina/alignment/masked_alignment.bcf.gz"
    shell:
        """
        bcftools view \
            --targets-file ^{input.mask} \
            --output-type b \
            --output {output.masked_alignment} \
            {input.alignment}
        """


def determine_alignment_from_vcf_input( wildcards ):
    if config["MASK"]:
        return rules.mask_vcf.output.masked_alignment
    return determine_mask_vcf_inputs( wildcards )


rule generate_alignment_from_vcf:
    input:
        vcf=determine_alignment_from_vcf_input,
        reference=rules.concatenate_reference.output.concatenated_reference
    output:
        fasta_alignment="intermediates/illumina/alignment/masked_alignment.fasta"
    shell:
        """
        custom-script \
            --vcf {input.vcf} \
            --reference {input.reference} \
            --output {output.fasta_alignment}
        """


rule run_gubbins:
    input:
        alignment=rules.generate_alignment_from_vcf.output.fasta_alignment
    params:
        prefix="intermediates/illumina/recombination_detection/gubbins",
        tree_builder="hybrid",
        substitution_model=config["tree_building"]["model"],
        gubbins_options=""
    threads: 16
    output:
        masked_alignment="intermediates/illumina/recombination_detection/gubbins.filtered_polymorphic_sites.fasta",
        recombinant_sites="intermediates/illumina/recombination_detection/gubbins.recombination_predictions.gff",
        other=temp( expand( "intermediates/illumina/recombination_detection/gubbins.{extension}",
            extension=[
                "recombination_predictions.embl",
                "branch_base_reconstruction.embl",
                "summary_of_snp_distribution.vcf",
                "per_branch_statistics.csv",
                "filtered_polymorphic_sites.phylip",
                "node_labelled.final_tree.tre",
                "log"
            ]
        ) )
    shell:
        """
        run_gubbins.py \
            --prefix {params.prefix} \
            --tree-builder {params.tree_builder} \
            --model {params.substitution_model} \
            --threads {threads} \
            {input.alignment}
        """


def calculate_outgroup( wildcards ):
    outgroup = config["tree_building"]["outgroup"]
    return f"-o {outgroup:q}" if outgroup != "" else ""


def determine_tree_input( wildcards ):
    if config["DETECT"]:
        return rules.run_gubbins.output.masked_alignment
    return rules.generate_alignment_from_vcf.output.fasta_alignment


rule generate_tree:
    input:
        alignment=determine_tree_input
    params:
        model=config["tree_building"]["model"],
        iqtree_parameters=config["tree_building"]["iqtree_parameters"],
        outgroup=calculate_outgroup
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
            -nt AUTO \
            -m {params.model} \
            {params.outgroup} \
            {params.iqtree_parameters} \
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
    #        reference=rules.concatenate_reference.output.concatenated_reference
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
