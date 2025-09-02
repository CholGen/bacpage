version 1.0

workflow reference_based_assembly {
    input {
        File read1
        File read2
        String sample_name
        File reference = "gs://bacpage-resources/vc_reference.fasta"

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

        Int disk_size = 8 # in GiB? Should check the size of the input.
        Int memory = 4
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
        File        consensus_sequence = ref_based_assembly.consensus_sequence
        File        samtools_idxstats = ref_based_assembly.samtools_idxstats
        File        samtools_stats = ref_based_assembly.samtools_stats
        File        fastqc_data = ref_based_assembly.fastqc_data
        File        bamqc_data = ref_based_assembly.bamqc_data
        File?       coverage_plot = ref_based_assembly.coverage_plot
        Float       total_reads = ref_based_assembly.total_reads
        Float       mapped_reads = ref_based_assembly.mapped_reads
        Float       percent_mapped_reads = ref_based_assembly.percent_mapped_reads
        Float       percent_coverage = ref_based_assembly.percent_coverage
        Int         median_depth = ref_based_assembly.median_depth
        String      bacpage_version = ref_based_assembly.bacpage_version
    }
}
task ref_based_assembly {
    input {
        File read1
        File read2
        String sample_name
        File reference

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

        Int disk_size = 8 # in GiB? Should check the size of the input.
        Int memory = 4
        Int cpu = 8
    }
    command <<<
        set -eux -o pipefail

        bacpage version > version

        bacpage setup tmp/
        cp ~{read1} tmp/input/~{sample_name}_R1.fastq.gz
        cp ~{read2} tmp/input/~{sample_name}_R2.fastq.gz

        # Generate the sample_data file so that weird names can be used.
        r1=$(realpath tmp/input/~{sample_name}_R1.fastq.gz)
        r2=$(realpath tmp/input/~{sample_name}_R2.fastq.gz)
        echo $"sample,read1,read2" > tmp/sample_data.csv
        echo "~{sample_name},${r1},${r2}" >> tmp/sample_data.csv

        cp ~{reference} reference.fasta

        ref=$(realpath reference.fasta)

        # TODO: generate the config.yaml from optional inputs.
        cat << EOF > tmp/config.yaml
        run_type: "Illumina"
        reference: $ref

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

        # Collect the stats!
        python << CODE
        import pandas as pd
        df = pd.read_csv( "tmp/results/reports/qc_report_data/multiqc_samtools_stats.txt", sep="\t" )
        with open( "total_reads", "w" ) as f: f.write( f"{df['sequences'][0]}\n" )
        with open( "mapped_reads", "w" ) as f: f.write( f"{df['reads_mapped'][0]}\n" )
        with open( "percent_mapped", "w" ) as f: f.write( f"{df['reads_mapped_percent'][0]:.1f}\n" )

        df = pd.read_csv( "tmp/results/reports/qc_report_data/multiqc_general_stats.txt", sep="\t" )
        print( df.columns )
        with open( "coverage", "w" ) as f: f.write( f"{df['qualimap_bamqc-10_x_pc'][0]:.1f}\n" )
        with open( "median_depth", "w" ) as f: f.write( f"{df['qualimap_bamqc-median_coverage'][0]:d}\n" )
        CODE

        # move results
        mv tmp/results/consensus/~{sample_name}.consensus.fasta ~{sample_name}.consensus.fasta
        mv tmp/results/reports/depth/~{sample_name}.depth.pdf ~{sample_name}.depth.pdf

        # save bamqc output to a tar.gz file
        tar --directory=tmp/results/reports/bamqc/ -czf ~{sample_name}_bamqc.tar.gz ~{sample_name}/

        # compress all qa_data
        # tar -czf ~{sample_name}.qa_data.tar.gz \
        #     ~{sample_name}_bamqc.tar.gz \
        #     tmp/results/reports/samtools/~{sample_name}.idxstats.txt \
        #     tmp/results/reports/samtools/~{sample_name}.stats.txt \
        #     tmp/results/reports/fastqc/~{sample_name}/~{sample_name}_fastqc.zip

    >>>
    output {
        File        consensus_sequence = "~{sample_name}.consensus.fasta"
        File        samtools_idxstats = "tmp/results/reports/samtools/~{sample_name}.idxstats.txt"
        File        samtools_stats = "tmp/results/reports/samtools/~{sample_name}.stats.txt"
        File        fastqc_data = "tmp/results/reports/fastqc/~{sample_name}/~{sample_name}_fastqc.zip"
        File        bamqc_data = "~{sample_name}_bamqc.tar.gz"
        File?       coverage_plot = "~{sample_name}.depth.pdf"
        Float       total_reads = read_float("total_reads")
        Float       mapped_reads = read_float( "mapped_reads" )
        Float       percent_mapped_reads = read_float( "percent_mapped" )
        Float       percent_coverage = read_float( "coverage" )
        Int         median_depth = read_int( "median_depth" )
        String      bacpage_version = read_string( "version" )

    }
    runtime {
        docker: "watronfire/bacpage:2025.08.21"
        cpu: cpu
        memory: memory + " GiB"
        disks: "local-disk " + disk_size + " HDD"
        bootDiskSizeGb: 50
    }
    meta {
        author: "Nathaniel L. Matteson"
        email: "nmatteson@bwh.harvard.edu"
        description: "## Reference-based assembly with bacpage \n This workflow performs referenced-based assembly using the bacpage pipeline."
    }
}


