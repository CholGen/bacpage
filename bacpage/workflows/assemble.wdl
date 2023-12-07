version 1.0

workflow reference_based_assembly {
    input {
        File read1
        File read2
        String sample_name
        File reference = "https://github.com/CholGen/bacpage/raw/split_into_command/bacpage/resources/vc_reference.fasta"

        String? bwa_parameters = "-M"

        Int? minimum_length = 30
        Int? minimum_quality = 20
        Int? window_length = 4

        Int? required_depth = 15

        Int? plotting_bin_size = 10000

        Int? maximum_depth = 2000
        Int? minimum_mapping_quality = 30
        Int? minimum_base_quality = minimum_quality
        String? mpileup_parameters = "-B -a INFO/AD,INFO/ADF,INFO/ADR -Ou"
        String? call_parameters = "-mv -Ov --ploidy 1"

        Int? minimum_depth = required_depth
        Float? minimum_support = 0.5
        Int? minimum_strand_depth = 1

        String? consensus_parameters = "--mark-del N"

        Int disk_size = 16 # in GiB? Should check the size of the input.
        Int memory = 16
        Int cpu = 8
    }
    call ref_based_assembly {
        input:
            read1 = read1,
            read2 = read2,
            sample_name = sample_name,
            reference = reference,
            bwa_parameters = bwa_parameters,
            minimum_length = minimum_length,
            minimum_quality = minimum_quality,
            window_length = window_length,
            required_depth = required_depth,
            plotting_bin_size = plotting_bin_size,
            maximum_depth = maximum_depth,
            minimum_mapping_quality = minimum_mapping_quality,
            minimum_base_quality = minimum_base_quality,
            mpileup_parameters = mpileup_parameters,
            call_parameters = call_parameters,
            minimum_depth = minimum_depth,
            minimum_support = minimum_support,
            minimum_strand_depth = minimum_strand_depth,
            consensus_parameters = consensus_parameters,
            disk_size = disk_size,
            memory = memory,
            cpu = cpu,
    }

    output {
        File consensus_sequence = ref_based_assembly.consensus_sequence
    }
}

task ref_based_assembly {
    input {
        File read1
        File read2
        String sample_name
        File reference = "vc_reference.fasta"

        String? bwa_parameters = "-M"

        Int? minimum_length = 30
        Int? minimum_quality = 20
        Int? window_length = 4

        Int? required_depth = 15

        Int? plotting_bin_size = 10000

        Int? maximum_depth = 2000
        Int? minimum_mapping_quality = 30
        Int? minimum_base_quality = minimum_quality
        String? mpileup_parameters = "-B -a INFO/AD,INFO/ADF,INFO/ADR -Ou"
        String? call_parameters = "-mv -Ov --ploidy 1"

        Int? minimum_depth = required_depth
        Float? minimum_support = 0.5
        Int? minimum_strand_depth = 1

        String? consensus_parameters = "--mark-del N"

        Int disk_size = 16 # in GiB? Should check the size of the input.
        Int memory = 16
        Int cpu = 8
    }
    command <<<
        bacpage setup tmp/
        echo "sample,read1,read2\n~{sample_name},~{read1},~{read2}" > tmp/sample_data.csv

        # TODO: generate the config.yaml from optional inputs.
        cat << EOF > tmp/config.yaml
        run_type: "Illumina"
        reference: ~{reference}

        preprocessing:
          check_size: False
          minimum_size: 100

        alignment_bwa:
          bwa_params: ~{bwa_parameters}

        trimming:
          minimum_length: ~{minimum_length}
          minimum_quality: ~{minimum_quality}
          window_length: ~{window_length}

        coverage_mask:
          required_depth: ~{required_depth}

        plot_coverage:
          bin_size: ~{plotting_bin_size}

        call_variants:
          maximum_depth: ~{maximum_depth}
          minimum_mapping_quality: ~{minimum_mapping_quality}
          minimum_base_quality: ~{minimum_base_quality}
          mpileup_parameters: ~{mpileup_parameters}
          call_parameters: ~{call_parameters}

        filter_variants:
          minimum_depth: ~{minimum_depth}
          minimum_support: ~{minimum_support}
          minimum_strand_depth: ~{minimum_strand_depth}

        call_consensus:
          consensus_parameters: ~{consensus_parameters}
        EOF

        bacpage assemble tmp/

        # zip results
        mv tmp/results/consensus/~{sample_name}.consensus.fasta consensus.fasta
    >>>
    output {
        File consensus_sequence = "tmp/results/consensus/~{sample_name}.consensus.fasta"
    }
    runtime {
        docker: "watronfire/bacpage:latest"
        cpu: cpu
        memory: memory + " GiB"
        disks: "local-disk " + disk_size + " HDD"
    }
}
