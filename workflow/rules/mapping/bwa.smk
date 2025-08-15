
rule bwa__build_index:
    input:
        "{reference_dir}/{fasta}.fa",
    output:
        idx=protected(multiext("{reference_dir}/bwa_index/{fasta}", ".amb", ".ann", ".bwt", ".pac", ".sa")),
    params:
        extra="-a bwtsw",
    wildcard_constraints:
        fasta="|".join(get_reference_names()),
    log:
        "{reference_dir}/bwa_index/logs/{fasta}.log",
    wrapper:
        "v7.2.0/bio/bwa/index"


rule bwa__map:
    input:
        reads=get_fastq_for_mapping,
        idx=infer_bwa_index_for_mapping,
    output:
        bam=temp_mapping("results/mapping/{reference}/{sample}.original.bam"),
        idx=temp_mapping("results/mapping/{reference}/{sample}.original.bam.bai"),
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:UNKNOWN'",
        sorting="samtools",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra="",  # Extra args for samtools/picard.
    threads: min(config["threads"]["mapping__mapping"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_mapping,
    log:
        "logs/mapping/bwa/{reference}/{sample}.log",
    wrapper:
        "v7.2.0/bio/bwa/mem"


ruleorder: bwa__map > samtools__bam_index
