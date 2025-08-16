rule quast__evaluate_assembly:
    input:
        fasta="results/assembly/{sample}/{assembly_tool}/contigs.fasta",
    output:
        pdf=report(
            "results/assembly/{sample}/{assembly_tool}/QUAST/report.pdf",
            category="Assembly",
            labels={
                "Sample": "{sample}",
                "Type": "QUAST report",
            },
        ),
        basic_reports=multiext("results/assembly/{sample}/{assembly_tool}/QUAST/report.", "html", "tex", "txt", "tsv"),
        transposed=multiext("results/assembly/{sample}/{assembly_tool}/QUAST/transposed_report.", "tex", "txt", "tsv"),
        basic_stats=directory("results/assembly/{sample}/{assembly_tool}/QUAST/basic_stats/"),
        icarus="results/assembly/{sample}/{assembly_tool}/QUAST/icarus.html",
        viewer="results/assembly/{sample}/{assembly_tool}/QUAST/icarus_viewers/contig_size_viewer.html",
        log="results/assembly/{sample}/{assembly_tool}/QUAST/quast.log",
    log:
        "logs/report/quast/{assembly_tool}/{sample}.log",
    threads: get_threads_for_report()
    params:
        extra=get_quast_params(),
    wrapper:
        "v7.2.0/bio/quast"
