rule bandage__visualise_contig_overlaps:
    input:
        "results/assembly/{sample}/{assembly_tool}/assembly_graph_with_scaffolds.gfa",
    output:
        report(
            "results/assembly/{sample}/{assembly_tool}/bandage/bandage.svg",
            category="Assembly",
            labels={
                "Sample": "{sample}",
                "Type": "Bandage SVG",
            },
        ),
    params:
        dir=lambda wildcards, output: os.path.dirname(output[0]),
    log:
        "logs/report/bandage_svg/{assembly_tool}/{sample}.log",
    conda:
        "../../envs/bandage.yaml"
    localrule: True
    shell:
        "(mkdir -p {params.dir} && Bandage image {input} {output}) > {log} 2>&1"


rule bandage__info:
    input:
        "results/assembly/{sample}/{assembly_tool}/assembly_graph_with_scaffolds.gfa",
    output:
        "results/assembly/{sample}/{assembly_tool}/bandage/bandage.info",
    params:
        dir=lambda wildcards, output: os.path.dirname(output[0]),
    log:
        "logs/report/bandage_info/{assembly_tool}/{sample}.log",
    conda:
        "../../envs/bandage.yaml"
    localrule: True
    shell:
        "(mkdir -p {params.dir} && Bandage info {input} > {output}) 2> {log}"
