version 1.0

workflow phylogeny_reconstruction {
    input {
        Array[File] consensus_sequences
        File reference = "https://github.com/CholGen/bacpage/raw/split_into_command/bacpage/resources/vc_reference.fastaclear"
        File? recombinant_mask = "https://github.com/CholGen/bacpage/raw/split_into_command/bacpage/resources/cholera_mask.gff"
        File? background_dataset = "https://github.com/watronfire/Alignment-Repo/raw/main/vibrio-cholerae/o1_cholera_v1.vcf.gz"

        Float? minimum_completeness = 0.9
        String? outgroup = ""
        String? model = "GTR"
        String? iqtree_parameters = "-nt AUTO -m TEST -bb 1000"

        Boolean skip_detection = false

        Int disk_size = 32 # in GiB? Should check the size of the input.
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
    }
}

task build_phylogeny {
    input {
        Array[File] consensus_sequences
        File reference = "https://github.com/CholGen/bacpage/raw/split_into_command/bacpage/resources/vc_reference.fasta"
        File? recombinant_mask = "https://github.com/CholGen/bacpage/raw/split_into_command/bacpage/resources/cholera_mask.gff"
        File? background_dataset = ""

        Float? minimum_completeness = 0.9
        String? outgroup = ""
        String? model = "GTR"
        String? iqtree_parameters = "-nt AUTO -m TEST -bb 1000"

        Boolean skip_detection = false

        Int disk_size = 32 # in GiB? Should check the size of the input.
        Int memory = 16
        Int cpu = 16
    }
    command <<<
        set -euxo pipefail

        mkdir tmp/
        cp ~{sep=" " consensus_sequences} tmp/

        # Construct config
        cat << EOF > tmp/config.yaml
        reference: ~{reference}
        recombinant_mask: "~{default="" recombinant_mask}"
        background_dataset: "~{default="" background_dataset}"

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
        if [ ~{!skip_detection} ] ; then
            mv tmp/results/phylogeny/recombinant_regions.gff .
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
