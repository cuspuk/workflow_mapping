rule multiqc__report:
    input:
        unpack(infer_multiqc_inputs_for_reference),
        config=f"{workflow.basedir}/resources/multiqc.yaml",
    output:
        report(
            "results/_aggregation/multiqc_{reference}.html",
            category="Summary",
            labels={
                "Reference": "{reference}",
                "Type": "MultiQC",
            },
        ),
    params:
        use_input_files_only=True,
        extra=lambda wildcards: f"--title 'Reference: {wildcards.reference}'",
    log:
        "logs/multiqc/{reference}.log",
    wrapper:
        "v3.13.8/bio/multiqc"
