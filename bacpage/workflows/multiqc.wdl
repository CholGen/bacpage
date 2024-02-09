version 1.0

workflow MultiQC {
    input {
        Array[File]?    samtools_idxstats
        Array[File]?    samtools_stats
        Array[File]?    fastqc_data
        Array[File]     bamqc_data
        Array[File]?    quast_reports
        Array[File]?    busco_reports

        Array[String]?  gambit_results
        Array[String]?  sample_ids

        Boolean         force = false
        Boolean         full_names = false
        String?         title
        String?         comment
        String?         file_name
        String          out_dir = "./multiqc-output"
        String?         template
        String?         tag
        String?         ignore_analysis_files
        String?         ignore_sample_names
        File?           sample_names
        Array[String]?  exclude_modules
        Array[String]?  module_to_use
        Boolean         data_dir = false
        Boolean         no_data_dir = false
        String?         output_data_format
        Boolean         zip_data_dir = false
        Boolean         export = false
        Boolean         flat = false
        Boolean         interactive = true
        Boolean         lint = false
        Boolean         pdf = false
        Boolean         megaQC_upload = false # Upload generated report to MegaQC if MegaQC options are found
        File?           config = "gs://bacpage-resources/multiqc_config.yaml"
        String?         config_yaml

        String          docker = "quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"
    }
    call MultiQC_task {
        input:
            samtools_idxstats = samtools_idxstats,
            samtools_stats = samtools_stats,
            fastqc_data = fastqc_data,
            bamqc_data = bamqc_data,
            quast_reports = quast_reports,
            busco_reports = busco_reports,
            gambit_results = gambit_results,
            sample_ids = sample_ids,
            force = force,
            full_names = full_names,
            title = title,
            comment = comment,
            file_name = file_name,
            out_dir = out_dir,
            template = template,
            tag = tag,
            ignore_analysis_files = ignore_analysis_files,
            ignore_sample_names = ignore_sample_names,
            sample_names = sample_names,
            exclude_modules = exclude_modules,
            module_to_use = module_to_use,
            data_dir = data_dir,
            no_data_dir = no_data_dir,
            output_data_format = output_data_format,
            zip_data_dir = zip_data_dir,
            export = export,
            flat = flat,
            interactive = interactive,
            lint = lint,
            pdf = pdf,
            megaQC_upload = megaQC_upload,
            config = config,
            config_yaml = config_yaml,
            docker = docker,
    }
    output {
      File multiqc_report = MultiQC_task.multiqc_report
      File multiqc_data_dir_tarball = MultiQC_task.multiqc_data_dir_tarball
      File multiqc_gambit_file = MultiQC_task.gambit_table
    }
}

task MultiQC_task {
    input {
        Array[File]?        samtools_idxstats
        Array[File]?        samtools_stats
        Array[File]?        fastqc_data
        Array[File]         bamqc_data
        Array[File]?        quast_reports
        Array[File]?        busco_reports

        Array[String]?      gambit_results
        Array[String]?      sample_ids

        Boolean             force = false
        Boolean             full_names = false
        String?             title
        String?             comment
        String?             file_name
        String              out_dir = "./multiqc-output"
        String?             template
        String?             tag
        String?             ignore_analysis_files
        String?             ignore_sample_names
        File?               sample_names
        Array[String]?      exclude_modules
        Array[String]?      module_to_use
        Boolean             data_dir = false
        Boolean             no_data_dir = false
        String?             output_data_format
        Boolean             zip_data_dir = false
        Boolean             export = false
        Boolean             flat = false
        Boolean             interactive = true
        Boolean             lint = false
        Boolean             pdf = false
        Boolean             megaQC_upload = false # Upload generated report to MegaQC if MegaQC options are found
        File?               config = "gs://bacpage-resources/multiqc_config.yaml"
        String?             config_yaml

        String              docker = "quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"
    }

    parameter_meta {
        output_data_format: { description: "[tsv|yaml|json] default:tsv" }
    }

    # get the basename in all wdl use the filename specified (sans ".html" extension, if specified)
    String report_filename = if (defined(file_name)) then basename(select_first([file_name]), ".html") else "multiqc"
    Int disk_size = 375

    command <<<
        set -ex -o pipefail

        # Move files
        mkdir tmp/

        cp ~{sep=" " bamqc_data} tmp/
        find tmp/ -name '*_bamqc.tar.gz' -execdir tar -xzf '{}' ';'

        if [ ~{defined(fastqc_data)} ]; then
            cp ~{sep=" " fastqc_data} tmp/
        fi

        if [ ~{defined(samtools_idxstats)} ]; then
            cp ~{sep=" " samtools_idxstats} tmp/
        fi

        if [ ~{defined(samtools_stats)} ]; then
            cp ~{sep=" " samtools_stats} tmp/
        fi

        if [ ~{defined(busco_reports)} ]; then
            cp ~{sep=" " busco_reports} tmp/
        fi

        if [ ~{defined(quast_reports)} ]; then
            cp ~{sep=" " quast_reports} tmp/
        fi

        if [ ~{defined(gambit_results) && ~{defined(sample_ids)} ]; then
            cat << EOF > tmp/gambit_mqc.tsv
            # plot_type: 'generalstats'
            # headers:
            #   gambit:
            #       title: "Species prediction"
            #       description: "Predicted taxonomic classification based on GAMBIT"
            Sample\tgambit
            EOF

            paste ~{write_lines(sample_ids)} ~{write_lines(gambit_results)} > tmp.tsv

            cat ~{gambit_file} >> tmp/gambit_mqc.tsv
        fi

        multiqc \
        --outdir "~{out_dir}" \
        ~{true="--force" false="" force} \
        ~{true="--fullnames" false="" full_names} \
        ~{"--title " + title} \
        ~{"--comment " + comment} \
        ~{"--filename " + file_name} \
        ~{"--template " + template} \
        ~{"--tag " + tag} \
        ~{"--ignore " + ignore_analysis_files} \
        ~{"--ignore-samples" + ignore_sample_names} \
        ~{"--sample-names " + sample_names} \
        ~{true="--exclude " false="" defined(exclude_modules)}~{sep=' --exclude ' select_first([exclude_modules,[]])} \
        ~{true="--module " false="" defined(module_to_use)}~{sep=' --module ' select_first([module_to_use,[]])} \
        ~{true="--data-dir" false="" data_dir} \
        ~{true="--no-data-dir" false="" no_data_dir} \
        ~{"--data-format " + output_data_format} \
        ~{true="--zip-data-dir" false="" zip_data_dir} \
        ~{true="--export" false="" export} \
        ~{true="--flat" false="" flat} \
        ~{true="--interactive" false="" interactive} \
        ~{true="--lint" false="" lint} \
        ~{true="--pdf" false="" pdf} \
        ~{false="--no-megaqc-upload" true="" megaQC_upload} \
        ~{"--config " + config} \
        ~{"--cl-config " + config_yaml } \
        tmp/

        if [ -z "~{file_name}" ]; then
            mv "~{out_dir}/~{report_filename}_report.html" "~{out_dir}/~{report_filename}.html"
        fi

        tar -c "~{out_dir}/~{report_filename}_data" | gzip -c > "~{report_filename}_data.tar.gz"
        >>>

    output {
        File multiqc_report           = "~{out_dir}/~{report_filename}.html"
        File multiqc_data_dir_tarball = "~{report_filename}_data.tar.gz"
        File gambit_table = "tmp/gambit_mqc.tsv"
    }

    runtime {
        memory: "8 GB"
        cpu: 16
        docker: "~{docker}"
        disks:  "local-disk " + disk_size + " LOCAL"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem2_ssd1_v2_x2"
        maxRetries: 2
    }
}
