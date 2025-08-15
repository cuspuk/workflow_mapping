
rule blast__download_database:
    output:
        blast_db=protected(
            multiext(
                "{blast_db_dir}/{reference_tag}.{type}",
                "db",
                "ot",
                "tf",
                "to",
            )
        ),
    params:
        blast_db_dir=lambda wildcards, output: os.path.dirname(output.blast_db[0]),
    log:
        "{blast_db_dir}/logs/{reference_tag}_{type}.log",
    retries: 3
    conda:
        "../../envs/blast.yaml"
    shell:
        "(mkdir -p {params.blast_db_dir} && cd {params.blast_db_dir} && update_blastdb.pl --decompress {wildcards.reference_tag} && blastdbcheck -db {wildcards.reference_tag} ) > {log} 2>&1"


rule blast__create_header:
    output:
        header=temp("results/classification/blast_header.tsv"),
    params:
        header="\t".join(BLAST_HEADER.split()),
    log:
        "logs/classification/blast/header.log",
    localrule: True
    conda:
        "../../envs/echo.yaml"
    shell:
        "echo -e {params.header:q}  > {output.header}"


rule blast__query:
    input:
        blast_db=get_blast_db(),
        contigs="results/assembly/{sample}/{assembly_tool}/contigs.fasta",
        header="results/classification/blast_header.tsv",
    output:
        tsv="results/classification/{sample}/{assembly_tool}/blast.tsv",
        tsv_headerless=temp("results/classification/{sample}/{assembly_tool}/blast.tsv.tmp"),
    params:
        blast_db_dir=lambda wildcards, input: os.path.dirname(input.blast_db[0]),
        binary=get_blast_binary(),
        blast_format="6 {header}".format(header=BLAST_HEADER),
        max_number_of_hits=get_max_number_of_hits(),
        reference_tag=get_blast_ref_tag(),
    threads: get_threads_for_classification()
    resources:
        mem_mb=get_mem_mb_for_classification,
    log:
        "logs/classification/blast/{assembly_tool}/{sample}/query.log",
    conda:
        "../../envs/blast.yaml"
    shell:
        "(export BLASTDB={params.blast_db_dir} && {params.binary} -db {params.reference_tag} -query {input.contigs} -out {output.tsv_headerless}"
        " -outfmt {params.blast_format:q} -num_threads {threads} -max_target_seqs {params.max_number_of_hits}"
        " && cat {input.header} {output.tsv_headerless} > {output.tsv}"
        " ) > {log} 2>&1"


rule custom__split_contigs:
    input:
        contigs="results/assembly/{sample}/{assembly_tool}/contigs.fasta",
    output:
        fasta_dir=directory("results/classification/{sample}/{assembly_tool}/blast/annotation/sequences/"),
        seqinfo="results/classification/{sample}/{assembly_tool}/blast/annotation/attributes/seqinfo.fa",
    params:
        fasta=lambda wildcards, output: os.path.join(output.fasta_dir, "all_contigs.fa"),
    log:
        "logs/classification/blast/{assembly_tool}/{sample}/split_contigs.log",
    localrule: True
    conda:
        "../../envs/blast_report.yaml"
    script:
        "../../scripts/summary_contigs.py"


rule custom__summary_blast:
    input:
        tsv="results/classification/{sample}/{assembly_tool}/blast.tsv",
    output:
        tsv="results/classification/{sample}/{assembly_tool}/blast/annotation/attributes/blast.tsv",
    params:
        min_query_coverage=config["assembly__classification__blast"]["min_query_coverage"],
        max_target_seqs=config["assembly__classification__blast"]["max_target_seqs"],
    log:
        "logs/classification/blast/{assembly_tool}/{sample}/summary.log",
    localrule: True
    conda:
        "../../envs/blast_report.yaml"
    script:
        "../../scripts/summary_blast.py"


rule custom__summary_html:
    input:
        template=os.path.join(workflow.basedir, "resources", "abundances.html.txt"),
        fasta_dir="results/classification/{sample}/{assembly_tool}/blast/annotation/sequences/",
        attrs_seqinfo="results/classification/{sample}/{assembly_tool}/blast/annotation/attributes/seqinfo.fa",
        attrs_blast="results/classification/{sample}/{assembly_tool}/blast/annotation/attributes/blast.tsv",
    output:
        html=report(
            "results/classification/{sample}/{assembly_tool}/blast/summary.html",
            category="Assembly",
            labels={
                "Sample": "{sample}",
                "Type": "Blast report",
            },
        ),
        tsv="results/classification/{sample}/{assembly_tool}/blast/summary.html.tsv",
    params:
        fasta=lambda wildcards, input: os.path.join(input.fasta_dir, "all_contigs.fa"),
        max_query_seqs=config["assembly__classification__blast"]["max_query_seqs"],
        seqs_per_page=config["assembly__classification__blast"]["seqs_per_page"],
        sort_by="Sequence",
        sort_how="asc",
        fasta_dir="sequences",
        html_columns=["Sequence", "Length", "Compress ratio", "Homologue link"],
    log:
        "logs/classification/blast/{assembly_tool}/{sample}/summary_html.log",
    localrule: True
    conda:
        "../../envs/blast_report.yaml"
    script:
        "../../scripts/summary.py"
