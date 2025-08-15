rule spades__assemble_reads_into_contigs:
    input:
        unpack(get_fastq_for_assembly),
    output:
        fastg="results/assembly/{sample}/spades/assembly_graph.fastg",
        gfa="results/assembly/{sample}/spades/assembly_graph_with_scaffolds.gfa",
        contigs="results/assembly/{sample}/spades/contigs.fasta",
        scaffolds="results/assembly/{sample}/spades/scaffolds.fasta",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.contigs),
        extra=get_spades_params(),
    threads: get_threads_for_assembly()
    resources:
        mem_mb=get_mem_mb_for_assembly,
    log:
        "logs/assembly/spades/{sample}.log",
    conda:
        "../../envs/spades.yaml"
    shell:
        "spades.py -1 {input.r1} -2 {input.r2} -o {params.outdir} --threads {threads} {params.extra} > {log} 2>&1"
