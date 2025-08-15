rule samtools__index_reference:
    input:
        reference="{reference_dir}/{reference}.fa",
    output:
        protected("{reference_dir}/{reference}.fa.fai"),
    log:
        "{reference_dir}/logs/samtools__index_reference/{reference}.log",
    params:
        extra="",
    wrapper:
        "v7.2.0/bio/samtools/faidx"


rule picard__prepare_dict_index:
    input:
        "{reference_dir}/{reference}.fa",
    output:
        protected("{reference_dir}/{reference}.dict"),
    log:
        "{reference_dir}/logs/picard__prepare_dict_index/{reference}.log",
    params:
        extra="",
    resources:
        mem_mb=get_mem_mb_for_deduplication,
    wrapper:
        "v7.2.0/bio/picard/createsequencedictionary"


rule samtools__bam_index:
    input:
        "results/mapping/{reference}/{sample}.original.bam",
    output:
        temp_mapping("results/mapping/{reference}/{sample}.original.bam.bai"),
    log:
        "logs/samtools_index/{reference}/original/{sample}.log",
    params:
        extra="",
    threads: min(config["threads"]["mapping__indexing"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_indexing,
    wrapper:
        "v7.2.0/bio/samtools/index"


rule samtools__stats:
    input:
        bam=infer_final_bam,
        ref=infer_reference_fasta,
    output:
        report(
            "results/mapping/{reference}/{sample}.{bam_step}.samtools_stats",
            category="Mapping QC for {reference}",
            labels={
                "Sample": "{sample}",
                "Type": "Samtools stats - {bam_step}",
            },
        ),
    params:
        extra=lambda wildcards, input: f"--ref-seq {input.ref}",
    threads: min(config["threads"]["mapping__indexing"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_indexing,
    log:
        "logs/mapping/samtools_stats/{reference}/{sample}_{bam_step}.log",
    wrapper:
        "v7.2.0/bio/samtools/stats"
