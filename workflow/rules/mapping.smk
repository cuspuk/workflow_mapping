rule samtools__index_reference:
    input:
        reference="{reference_dir}/{reference}.fa",
    output:
        protected("{reference_dir}/{reference}.fa.fai"),
    log:
        "{reference_dir}/logs/samtools__prepare_fai_index/{reference}.log",
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.13.4/wrappers/samtools/faidx"


rule custom__infer_read_group:
    input:
        get_fastq_for_mapping,
    output:
        read_group="results/reads/.read_groups/{sample}.txt",
    params:
        sample_id=lambda wildcards: wildcards.sample,
    log:
        "logs/custom/infer_and_store_read_group/{sample}.log",
    localrule: True
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.13.4/wrappers/custom/read_group"


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
        "v3.13.7/bio/picard/createsequencedictionary"


rule samtools__bam_index:
    input:
        bam="results/mapping/{reference}/{sample}.original.bam",
    output:
        bai=temp_mapping("results/mapping/{reference}/{sample}.original.bam.bai"),
    threads: min(config["threads"]["mapping__indexing"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_indexing,
    log:
        "logs/mapping/indexing/{reference}/mapped/{sample}.log",
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.13.4/wrappers/samtools/index"


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
        "v3.13.7/bio/samtools/stats"
