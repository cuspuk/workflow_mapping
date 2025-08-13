rule minimap2_index:
    input:
        target="{reference_dir}/{fasta}.fa",
    output:
        protected("{reference_dir}/{fasta}.mmi"),
    log:
        "{reference_dir}/minimap2_index/logs/{fasta}.log",
    params:
        extra="",
    threads: min(config["threads"]["mapping__mapping"], config["max_threads"])
    wildcard_constraints:
        fasta="|".join(get_reference_names()),
    wrapper:
        "v7.2.0/bio/minimap2/index"


rule minimap2_bam_sorted:
    input:
        target=infer_minimap2_index_for_mapping,
        query=get_fastq_for_mapping,
    output:
        bam=temp_mapping("results/mapping/{reference}/{sample}.original.bam"),
    log:
        "logs/mapping/minimap2/{reference}/{sample}.log",
    params:
        extra=lambda wildcards: f"{config['mapping__mapping__minimap2']['preset']} -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}'",
        sorting="coordinate",  # 'none', 'queryname' or 'coordinate'
        sort_extra="",  # optional: extra arguments for samtools/picard
    threads: min(config["threads"]["mapping__mapping"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_mapping,
    wrapper:
        "v7.2.0/bio/minimap2/aligner"
