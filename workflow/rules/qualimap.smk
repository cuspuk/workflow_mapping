rule qualimap__mapping_quality_report:
    input:
        bam="results/mapping/{reference}/{sample}.{bam_step}.bam",
        bai="results/mapping/{reference}/{sample}.{bam_step}.bam.bai",
    output:
        report_dir=report(
            directory("results/mapping/{reference}/bamqc/{sample}.{bam_step}"),
            category="Mapping QC for {reference}",
            labels={
                "Sample": "{sample}",
                "Type": "Qualimap for {bam_step}",
            },
            htmlindex="qualimapReport.html",
        ),
    params:
        extra=[
            "--paint-chromosome-limits",
            "-outformat PDF:HTML",
        ],
    resources:
        mem_mb=get_mem_mb_for_qualimap,
    log:
        "logs/qualimap/{reference}/{bam_step}/{sample}.log",
    wrapper:
        "v7.2.0/bio/qualimap/bamqc"
