rule minion_qc:
    input:
        summary = "sequencing_summary.txt"
    params:
        image_format = config["minion_qc"]["format"]
    output:
        report = directory( "results/reports/minion_qc/" )
    shell:
        """
        Rscript MinIONQC.R \
            --input={input.summary} \
            --outputdirectory={output.report} \
            --format={params.image_format}
        """

