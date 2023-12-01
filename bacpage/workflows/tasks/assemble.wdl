version 1.0

task ref_based_assembly_batch {
    input {
        File sample_data # Need to check if I can supply multiple samples.
        File config

        Int disk_size = 25 # in GiB? Should check the size of the input.
        Int memory = 16
    }
    command <<<
        # TODO will have to automatically generate the config file.

        bacpage assemble --samples ~{sample_data} --configfile ~{config} .

        # zip qc_report
        zip -r qc_report.zip results/report/qc_report/

        # zip results
        zip -r sequences.zip results/consensus/*.consensus.fasta
    >>>
    output {
        File sequences = "sequences.zip"
        File report = "qc_report.html"
        File report_data = "qc_report.zip"
    }
    runtime {
        docker: "watronfire/bacpage:latest"
        cpu: "16"
        memory: memory + " GiB"
        disks: "local-disk " + disk_size + " HDD"
    }
}

task ref_based_assembly {
    input {
        File read1
        File read2
        File config
        String sample_name

        Int disk_size = 16 # in GiB? Should check the size of the input.
        Int memory = 16
        Int cpu = 8
    }
    command <<<
        bacpage setup tmp/
        echo "sample,read1,read2\n~{sample_name},~{read1},~{read2}" > tmp/sample_data.csv

        # TODO: generate the config.yaml from optional inputs.

        bacpage assemble --configfile ~{config} tmp/

        # zip results
        mv results/consensus/~{sample_name}.consensus.fasta consensus.fasta
    >>>
    output {
        File consensus_sequence = "consensus.fasta"
    }
    runtime {
        docker: "watronfire/bacpage:latest"
        cpu: cpu
        memory: memory + " GiB"
        disks: "local-disk " + disk_size + " HDD"
    }
}
