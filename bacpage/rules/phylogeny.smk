def determine_output( wildcards ):
    if config["TERRA"]:
        return determine_alignment_from_vcf_input( wildcards )
    outputs = ["results/phylogeny/phylogeny.tree", "results/phylogeny/sparse_alignment.fasta"]
    if config["DETECT"] and not isinstance( config["DETECT"],str ):
        outputs.append( "results/phylogeny/recombinant_regions.gff" )
    return outputs


rule all:
    input:
        determine_output


def calculate_complete_sequences( wildcards ):
    complete_sequences = list( config["SAMPLES"].values() )
    if config["BACKGROUND"] != "":
        if not config["BACKGROUND"].name.endswith( (".vcf", ".vcf.gz", ".bcf", ".bcf.gz") ):
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
        SEDOPTION=
        if [[ "$OSTYPE" == "darwin"* ]]; then
          SEDOPTION="-i ''"
        fi
        
        sed '1h;/>/d;H;$!d;x;s/\\n/@/;s/\\n//g;s/@/\\n/' {input.reference} |\
        sed $SEDOPTION -e '$a\\' > {output.concatenated_reference}
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
        REFERENCE=$(head -n1 {input.reference} | cut -f2 -d \> | cut -f1 -d" ") &&\
        echo "1 ${{REFERENCE}}" > {output.chromosome_name} &&\
        cat {input.reference} {input.alignment} > {output.temp_alignment} &&\
        snp-sites -v {output.temp_alignment} | bcftools annotate --samples ^'${{REFERENCE}}' --rename-chrs {output.chromosome_name} -O b -o {output.vcf} &&\
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
            --missing-to-ref \
            {input.background_sequences} {input.user_sequences} &&\
        bcftools index {output.combined_vcf}
        """


rule convert_gff_to_bed:
    message: "Coverts GFF file to BED file readable by BCFtools."
    input:
        mask=config["MASK_LOCATION"],
        reference=rules.concatenate_reference.output.concatenated_reference
    output:
        mask_bed="intermediates/illumina/alignment_mask.txt"
    run:
        import pandas as pd

        reference_name = shell( "head -n1 {input.reference} | cut -f2 -d \> | cut -f1 -d' '",read=True ).strip()
        df = pd.read_csv( input.mask,sep="\t",header=None )
        df["CHROM"] = reference_name
        df = df[["CHROM", 3, 4]]
        df.columns = ["CHROM", "BEG", "END"]
        df.to_csv( output.mask_bed,header=False,index=False,sep="\t" )


def determine_mask_vcf_inputs( wildcards ):
    if (config["BACKGROUND"] == "") or (config["BACKGROUND"].suffix in [".fa", ".fasta"]):
        return rules.convert_to_vcf.output.vcf
    return rules.combine_sequences_and_background_vcf.output.combined_vcf


rule mask_vcf:
    input:
        alignment=determine_mask_vcf_inputs,
        mask=rules.convert_gff_to_bed.output.mask_bed
    output:
        masked_alignment="intermediates/illumina/alignment/masked_alignment.bcf.gz"
    shell:
        """
        bcftools view \
            --targets-file ^{input.mask} \
            --output-type b \
            --output {output.masked_alignment} \
            {input.alignment} &&\
        bcftools index {output.masked_alignment}
        """


def determine_alignment_from_vcf_input( wildcards ):
    if config["MASK"]:
        return rules.mask_vcf.output.masked_alignment
    return determine_mask_vcf_inputs( wildcards )


rule generate_alignment_from_vcf:
    input:
        vcf=determine_alignment_from_vcf_input,
        reference=rules.concatenate_reference.output.concatenated_reference
    params:
        script_location=workflow.source_path( "../scripts/vcf_to_fasta.py" )
    output:
        fasta_alignment=temp( "intermediates/illumina/alignment/masked_alignment.fasta" )
    shell:
        """
        python {params.script_location} \
            --vcf {input.vcf} \
            --reference {input.reference} \
            --output {output.fasta_alignment}
        """


rule bypass_gubbins:
    input:
        vcf=determine_alignment_from_vcf_input,
        mask=lambda wildcards: config["DETECT"] if isinstance( config["DETECT"],str ) else ""
    output:
        masked_vcf="intermediates/illumina/recombination_detection/gubbins_masked.bcf.gz"
    shell:
        """
        bcftools view \
            --targets-file ^{input.mask} \
            --output-type b \
            --output {output.masked_vcf} \
            {input.vcf}
        bcftools index {output.masked_vcf}
        """


rule generate_alignment_from_bypassed_gubbins:
    input:
        masked_vcf=rules.bypass_gubbins.output.masked_vcf,
        reference=rules.concatenate_reference.output.concatenated_reference
    output:
        fasta_alignment=temp( "intermediates/illumina/recombination_detection/gubbins_masked_skipped.fasta" )
    params:
        script_location=workflow.source_path( "../scripts/vcf_to_fasta.py" )
    shell:
        """
        python {params.script_location} \
            --vcf {input.masked_vcf} \
            --reference {input.reference} \
            --output {output.fasta_alignment}
        """


rule run_gubbins:
    input:
        alignment=rules.generate_alignment_from_vcf.output.fasta_alignment
    params:
        prefix="intermediates/illumina/recombination_detection/gubbins",
        tree_builder="hybrid",
        substitution_model=config["tree_building"]["model"] + "GAMMA",
        gubbins_options=""
    threads: 16
    output:
        masked_alignment=temp( "intermediates/illumina/recombination_detection/gubbins.filtered_polymorphic_sites.fasta" ),
        recombinant_sites="intermediates/illumina/recombination_detection/gubbins.recombination_predictions.gff",
        other=temp(
            expand(
                "intermediates/illumina/recombination_detection/gubbins.{extension}",
                extension=[
                    "recombination_predictions.embl",
                    "branch_base_reconstruction.embl",
                    "summary_of_snp_distribution.vcf",
                    "per_branch_statistics.csv",
                    "filtered_polymorphic_sites.phylip",
                    "node_labelled.final_tree.tre",
                    "log"
                ]
            )
        )
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
        if isinstance( config["DETECT"],str ):
            return rules.generate_alignment_from_bypassed_gubbins.output.fasta_alignment
        return rules.run_gubbins.output.masked_alignment
    return rules.generate_alignment_from_vcf.output.fasta_alignment


rule sparsify_alignment:
    input:
        alignment=determine_tree_input
    output:
        alignment="results/phylogeny/sparse_alignment.fasta"
    shell:
        """
        snp-sites -o {output.alignment} {input.alignment}
        """


rule generate_tree:
    input:
        alignment=rules.sparsify_alignment.output.alignment
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
        cp {input.iqtree_output} {output.final_tree}
        """


rule move_recombinant_mask:
    input:
        mask=rules.run_gubbins.output.recombinant_sites
    output:
        mask="results/phylogeny/recombinant_regions.gff"
    shell:
        """
        cp {input.mask} {output.mask}
        """
