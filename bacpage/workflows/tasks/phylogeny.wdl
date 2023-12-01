version 1.0

task phylogeny_reconstruction {
    input {
        Array[File] consensus_sequences
        File config

        String? background  # Unused at the moment. Will condition bacpage options
        String? mask        # Unused at the moment. Will condition bacpage options

        Float? minimum_completeness = 0.9
        Boolean? skip_detection = false

        Int disk_size = 32 # in GiB? Should check the size of the input.
        Int memory = 16
        Int cpu = 16
    }
    command <<<
        # TODO Definitely need to do some setup here.
        mkdir project_directory/
        cp ~{consensus_sequences} project_directory/

        bacpage phylogeny --configfile ~{config} --minimum-completeness ~{minimum_completeness} project_directory/

        # Move output
        mv project_directory/results/phylogeny/phylogeny.tree project_directory/results/phylogeny/sparse_alignment.fasta project_directory/results/phylogeny/recombinant_regions.gff .
    >>>
    output {
        File phylogeny = "phylogeny.tree"
        File sparse_alignment = "sparse_alignment.fasta"
        File recombinant_regions = "recombinant_regions.gff"
    }
    runtime {
        docker: "watronfire/bacpage:latest"
        cpu: cpu
        memory: memory + " GiB"
        disks: "local-disk " + disk_size + " HDD"
    }
}
