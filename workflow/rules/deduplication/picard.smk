rule picard__mark_duplicates:
    input:
        bams="results/mapping/{reference}/{sample}.original.bam",
        bai="results/mapping/{reference}/{sample}.original.bam.bai",
    output:
        bam="results/mapping/{reference}/{sample}.deduplication.bam",
        idx="results/mapping/{reference}/{sample}.deduplication.bam.bai",
        metrics="results/mapping/{reference}/{sample}.deduplication.stats",
    resources:
        mem_mb=get_mem_mb_for_deduplication,
    log:
        "logs/deduplication/picard/{reference}/{sample}.log",
    wrapper:
        "v4.3.0/bio/picard/markduplicates"
