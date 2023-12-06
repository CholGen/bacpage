version 1.0

task phylogeny_reconstruction {
    input {
        Array[File] consensus_sequences
        File reference = "https://github.com/CholGen/bacpage/raw/split_into_command/bacpage/resources/vc_reference.fasta"
        File recombinant_mask = "https://github.com/CholGen/bacpage/raw/split_into_command/bacpage/resources/cholera_mask.gff"
        File background_dataset = ""

        Float? minimum_completeness = 0.9
        String? outgroup = ""
        String? model = "GTR"
        String? iqtree_parameters = "-nt AUTO -m TEST -bb 1000"

        Int disk_size = 32 # in GiB? Should check the size of the input.
        Int memory = 16
        Int cpu = 16
    }
    command <<<
        # TODO Definitely need to do some setup here.
        mkdir tmp/
        cp ~{consensus_sequences} tmp/

        # Construct config
        cat << EOF > tmp/config.yaml
        reference: ~{reference}
        recombinant_mask: ~{recombinant_mask}
        background_dataset: ~{background_dataset}

        tree_building:
          minimum_completeness: ~{minimum_completeness}
          outgroup: ~{outgroup}
          model: ~{model}
          iqtree_parameters: ~{iqtree_parameters}
        EOF

        bacpage phylogeny project_directory/

        # Move output
        mv project_directory/results/phylogeny/phylogeny.tree project_directory/results/phylogeny/sparse_alignment.fasta project_directory/results/phylogeny/recombinant_regions.gff .
    >>>
    output {
        File phylogeny = "phylogeny.tree"
        File sparse_alignment = "sparse_alignment.fasta"
        File? recombinant_regions = "recombinant_regions.gff"
    }
    runtime {
        docker: "watronfire/bacpage:latest"
        cpu: cpu
        memory: memory + " GiB"
        disks: "local-disk " + disk_size + " HDD"
    }
}
