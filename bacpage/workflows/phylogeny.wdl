version 1.0

workflow phylogeny_reconstruction {
    input {
        Array[File] consensus_sequences
        File background_dataset
        File reference = "gs://andersen-lab_temp/vc_reference.fasta"
        File? recombinant_mask

        Float? minimum_completeness = 0.9
        String? outgroup = ""
        String? model = "GTR"
        String? iqtree_parameters = "-nt AUTO -m TEST -bb 1000"

        Boolean skip_detection = false

        Int disk_size = 64 # in GiB? Should check the size of the input.
        Int memory = 16
        Int cpu = 16
    }
    call build_phylogeny {
        input:
            consensus_sequences = consensus_sequences,
            reference = reference,
            recombinant_mask = recombinant_mask,
            background_dataset = background_dataset,
            minimum_completeness = minimum_completeness,
            outgroup = outgroup,
            model = model,
            iqtree_parameters = iqtree_parameters,
            skip_detection=skip_detection,
            disk_size = disk_size,
            memory = memory,
            cpu = cpu,
    }
    output {
        File phylogeny = build_phylogeny.phylogeny
        File sparse_alignment = build_phylogeny.sparse_alignment
        File? recombinant_regions = build_phylogeny.recombinant_regions
        File? alignment = "combined_alignment.bcf.gz"
    }
}

task build_phylogeny {
    input {
        Array[File] consensus_sequences
        File background_dataset
        File reference
        File? recombinant_mask

        Float? minimum_completeness = 0.9
        String? outgroup = ""
        String? model = "GTR"
        String? iqtree_parameters = "-bb 1000"

        Boolean skip_detection = false

        Int disk_size = 64 # in GiB? Should check the size of the input.
        Int memory = 16
        Int cpu = 16
    }

    String background_name = basename( background_dataset )

    command <<<
        set -euxo pipefail

        mkdir tmp/
        cp ~{sep=" " consensus_sequences} tmp/

        cp ~{background_dataset} .
        background="$(pwd)/~{background_name}"

        # Construct config
        cat << EOF > tmp/config.yaml
        reference: ~{reference}
        recombinant_mask: "~{default="" recombinant_mask}"
        background_dataset: "$background"

        tree_building:
          minimum_completeness: ~{minimum_completeness}
          outgroup: "~{default="" outgroup}"
          model: ~{model}
          iqtree_parameters: ~{iqtree_parameters}
        EOF

        bacpage phylogeny \
            ~{true='--no-detect' false='' skip_detection} \
            tmp/

        # Move output
        mv tmp/results/phylogeny/phylogeny.tree tmp/results/phylogeny/sparse_alignment.fasta .
        mv tmp/intermediates/illumina/alignment/combined_alignment.bcf.gz .
        if [ ~{!skip_detection} ] ; then
            mv tmp/results/phylogeny/recombinant_regions.gff .
        fi
    >>>
    output {
        File phylogeny = "phylogeny.tree"
        File sparse_alignment = "sparse_alignment.fasta"
        File? recombinant_regions = "recombinant_regions.gff"
        File? alignment = "combined_alignment.bcf.gz"
    }
    runtime {
        docker: "watronfire/bacpage:latest"
        cpu: cpu
        memory: memory + " GiB"
        disks: "local-disk " + disk_size + " HDD"
    }
    meta {
        author: "Nathaniel L. Matteson"
        email: "nmatteson@bwh.harvard.edu"
        description: "## Phylogenetic Reconstruction with bacpage \n This workflow performs a phylogenetic reconstruction using the bacpage pipeline."
    }
}
