rule hisat2_indexL:
    input:
        fasta="{reference_dir}/{fasta}.fa",
    output:
        protected(
            multiext(
                "{reference_dir}/hisat2_index/{fasta}",
                ".1.ht2l",
                ".2.ht2l",
                ".3.ht2l",
                ".4.ht2l",
                ".5.ht2l",
                ".6.ht2l",
                ".7.ht2l",
                ".8.ht2l",
            )
        ),
    params:
        extra="--large-index",
    wildcard_constraints:
        fasta="|".join(get_reference_names()),
    log:
        "{reference_dir}/hisat2_indexL/logs/{fasta}.log",
    threads: min(config["threads"]["mapping__mapping"], config["max_threads"])
    wrapper:
        "v7.2.0/bio/hisat2/index"


rule hisat2_alignL:
    input:
        reads=get_fastq_for_mapping,
        idx=infer_hisat2_index_for_mapping,
    output:
        "results/mapping/{reference}/{sample}.unsorted.bam",
    log:
        "logs/mapping/hisat2/{reference}/{sample}.log",
    params:
        extra="",
    threads: min(config["threads"]["mapping__mapping"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_mapping,
    wrapper:
        "v7.2.0/bio/hisat2/align"


rule hisat2_samtools_sort:
    input:
        "results/mapping/{reference}/{sample}.unsorted.bam",
    output:
        temp_mapping("results/mapping/{reference}/{sample}.original.bam"),
    log:
        "logs/mapping/hisat2_sort/{reference}/{sample}.log",
    params:
        extra="",
    threads: min(config["threads"]["mapping__mapping"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_mapping,
    wrapper:
        "v7.2.0/bio/samtools/sort"
