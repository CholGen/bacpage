task MultiQC {
  input {
    Array[File]    input_files

    Boolean        force = false
    Boolean        full_names = false
    String?        title
    String?        comment
    String?        file_name
    String         out_dir = "./multiqc-output"
    String?        template
    String?        tag
    String?        ignore_analysis_files
    String?        ignore_sample_names
    File?          sample_names
    Array[String]? exclude_modules
    Array[String]? module_to_use
    Boolean        data_dir = false
    Boolean        no_data_dir = false
    String?        output_data_format
    Boolean        zip_data_dir = false
    Boolean        export = false
    Boolean        flat = false
    Boolean        interactive = true
    Boolean        lint = false
    Boolean        pdf = false
    Boolean        megaQC_upload = false # Upload generated report to MegaQC if MegaQC options are found
    File?          config = "gs://bacpage-resources/multiqc_config.yaml"
    String?        config_yaml

    String         docker = "quay.io/biocontainers/multiqc:1.8--py_2"
  }

  parameter_meta {
    output_data_format: { description: "[tsv|yaml|json] default:tsv" }
  }

  # get the basename in all wdl use the filename specified (sans ".html" extension, if specified)
  String report_filename = if (defined(file_name)) then basename(select_first([file_name]), ".html") else "multiqc"
  Int disk_size = 375

  command {
      set -ex -o pipefail

      echo "${sep='\n' input_files}" > input-filenames.txt
      echo "" >> input-filenames.txt

      multiqc \
      --file-list input-filenames.txt \
      --outdir "${out_dir}" \
      ${true="--force" false="" force} \
      ${true="--fullnames" false="" full_names} \
      ${"--title " + title} \
      ${"--comment " + comment} \
      ${"--filename " + file_name} \
      ${"--template " + template} \
      ${"--tag " + tag} \
      ${"--ignore " + ignore_analysis_files} \
      ${"--ignore-samples" + ignore_sample_names} \
      ${"--sample-names " + sample_names} \
      ${true="--exclude " false="" defined(exclude_modules)}${sep=' --exclude ' select_first([exclude_modules,[]])} \
      ${true="--module " false="" defined(module_to_use)}${sep=' --module ' select_first([module_to_use,[]])} \
      ${true="--data-dir" false="" data_dir} \
      ${true="--no-data-dir" false="" no_data_dir} \
      ${"--data-format " + output_data_format} \
      ${true="--zip-data-dir" false="" zip_data_dir} \
      ${true="--export" false="" export} \
      ${true="--flat" false="" flat} \
      ${true="--interactive" false="" interactive} \
      ${true="--lint" false="" lint} \
      ${true="--pdf" false="" pdf} \
      ${false="--no-megaqc-upload" true="" megaQC_upload} \
      ${"--config " + config} \
      ${"--cl-config " + config_yaml }

      if [ -z "${file_name}" ]; then
        mv "${out_dir}/${report_filename}_report.html" "${out_dir}/${report_filename}.html"
      fi

      tar -c "${out_dir}/${report_filename}_data" | gzip -c > "${report_filename}_data.tar.gz"
  }

  output {
      File multiqc_report           = "${out_dir}/${report_filename}.html"
      File multiqc_data_dir_tarball = "${report_filename}_data.tar.gz"
  }

  runtime {
    memory: "8 GB"
    cpu: 16
    docker: "${docker}"
    disks:  "local-disk " + disk_size + " LOCAL"
    disk: disk_size + " GB" # TES
    dx_instance_type: "mem2_ssd1_v2_x2"
    maxRetries: 2
  }
}
