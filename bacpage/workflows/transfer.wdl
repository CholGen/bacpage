version 1.0

workflow transfer_pe_reads {
    input {
        Array[String] read1
        Array[String] read1
        String target_bucket
    }
    call transfer_files as transfer_read1 {
        input:
            files_to_transfer = read1,
            target_bucket = target_bucket
    }
    call transfer_files as transfer_read2 {
        input:
            files_to_transfer = read1,
            target_bucket = target_bucket
    }
    output {
        File transferred_read1 = transfer_read1.transferred_files
        File transferred_read2 = transfer_read2.transferred_files
    }
}

task transfer_files {
    input {
        Array[String] files_to_transfer
        String target_bucket
        Int cpu = 4
        Int memory = 8
        String docker_image = "us-docker.pkg.dev/general-theiagen/theiagen/utility:1.1"
        Int disk_size = 100
    }
    meta {
        # added so that call caching is always turned off
        volatile: true
    }
    command <<<
        file_path_array="~{sep=' ' files_to_transfer}"

        gsutil -m cp -n ${file_path_array[@]} ~{target_bucket}

        echo "transferred_files" > transferred_files.tsv
        gsutil -m ls ~{target_bucket} >> transferred_files.tsv
     >>>
    output {
        File transferred_files = "transferred_files.tsv"
    }
    runtime {
        docker: "~{docker_image}"
        memory: "~{memory} GB"
        cpu: cpu
        disks: "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB"
        preemptible: 0
    }
}
