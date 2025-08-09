rule samtools_sort_queryname:
    input:
        bam="results/mapping/{reference}/{sample}.original.bam",
        bai="results/mapping/{reference}/{sample}.original.bam.bai",
    output:
        temp("results/mapping/{reference}/{sample}.queryname_sorted.bam"),
    log:
        "logs/deduplication/samtools_sort/{reference}/{sample}.log",
    threads: min(config["threads"]["mapping__mapping"], config["max_threads"])
    params:
        extra="-n",
    resources:
        mem_mb=get_mem_mb_for_deduplication,
    wrapper:
        "v7.2.0/bio/samtools/sort"


rule samtools_fixmate:
    input:
        bam="results/mapping/{reference}/{sample}.queryname_sorted.bam",
    output:
        bam=temp("results/mapping/{reference}/{sample}.fixmate.bam"),
    log:
        "logs/deduplication/samtools_fixmate/{reference}/{sample}.log",
    threads: min(config["threads"]["mapping__mapping"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_deduplication,
    params:
        extra="-m",
    wrapper:
        "v7.2.0/bio/samtools/fixmate/"


rule samtools_sort_after_fixmate:
    input:
        "results/mapping/{reference}/{sample}.fixmate.bam",
    output:
        temp("results/mapping/{reference}/{sample}.fixmate_sorted.bam"),
    log:
        "logs/deduplication/samtools_sort/{reference}/{sample}.log",
    threads: min(config["threads"]["mapping__mapping"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_deduplication,
    wrapper:
        "v7.2.0/bio/samtools/sort"


rule samtools_markdup:
    input:
        "results/mapping/{reference}/{sample}.fixmate_sorted.bam",
    output:
        bam="results/mapping/{reference}/{sample}.deduplication.bam",
        metrics="results/mapping/{reference}/{sample}.deduplication.stats",
    log:
        "logs/deduplication/samtools_markdup/{reference}/{sample}.log",
    params:
        extra="-c --no-PG",
    threads: min(config["threads"]["mapping__mapping"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_deduplication,
    wrapper:
        "v7.2.0/bio/samtools/markdup"


rule samtools_index:
    input:
        "results/mapping/{reference}/{sample}.deduplication.bam",
    output:
        "results/mapping/{reference}/{sample}.deduplication.bam.bai",
    log:
        "logs/deduplication/samtools_index/{reference}/{sample}.log",
    params:
        extra="",
    threads: min(config["threads"]["mapping__indexing"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_indexing,
    wrapper:
        "v7.2.0/bio/samtools/index"
