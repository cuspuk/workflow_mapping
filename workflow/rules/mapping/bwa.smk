
rule bwa__build_index:
    input:
        "{reference_dir}/{fasta}.fa",
    output:
        idx=protected(multiext("{reference_dir}/bwa_index/{fasta}", ".amb", ".ann", ".bwt", ".pac", ".sa")),
    params:
        prefix=lambda wildcards, output: os.path.splitext(output.idx[0])[0],
        approach="bwtsw",
    wildcard_constraints:
        fasta="|".join(get_reference_names()),
    log:
        "{reference_dir}/bwa_index/logs/{fasta}.log",
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.13.4/wrappers/bwa/index"


rule bwa__map_reads:
    input:
        reads=get_fastq_for_mapping,
        index=infer_bwa_index_for_mapping,
        read_group="results/reads/.read_groups/{sample}.txt",
    output:
        bam=temp_mapping("results/mapping/{reference}/{sample}.original.bam"),
    threads: min(config["threads"]["mapping__mapping"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_mapping,
    log:
        "logs/mapping/bwa/{reference}/{sample}.log",
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.13.4/wrappers/bwa/map"
