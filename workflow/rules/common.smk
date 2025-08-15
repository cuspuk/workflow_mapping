from snakemake.utils import validate


configfile: "config/config.yaml"


validate(config, "../schemas/config.schema.yaml", set_default=False)


### Layer for adapting other workflows  ###############################################################################


def get_fastq_for_mapping(wildcards):
    if config["mapping"]["_input"] == "reads":
        return reads_workflow.get_final_fastq_for_sample(wildcards.sample)
    elif config["mapping"]["_input"] == "assembly":
        return get_assembly_fasta(wildcards.sample)
    else:
        raise ValueError(f"Invalid input type for mapping: {config['mapping']['_input']}")


def get_fastq_for_assembly(wildcards):
    reads = reads_workflow.get_final_fastq_for_sample(wildcards.sample)
    return {
        "r1": reads[0],
        "r2": reads[1],
    }


def get_sample_names():
    return reads_workflow.get_sample_names()


def get_multiqc_inputs_for_reads():
    return reads_workflow.get_multiqc_inputs()


def get_reads_outputs():
    return reads_workflow.get_outputs()


### Data input handling independent of wildcards ######################################################################


def get_reference_name_from_path(path: str):
    name, ext = os.path.splitext(os.path.basename(os.path.realpath(path)))
    if ext not in [".fasta", ".fa"]:
        raise ValueError(f"Reference file {path} does not have a valid extension (.fasta or .fa)")
    return name


reference_dict = {
    get_reference_name_from_path(path): os.path.realpath(path) for path in config["mapping"]["reference_fasta_paths"]
}


def get_reference_dir_for_name(name: str):
    return os.path.dirname(reference_dict[name])


def get_reference_names():
    return reference_dict.keys()


def parse_blast_tag(path: str):
    return os.path.basename(os.path.dirname(os.path.realpath(path)))


def validate_blast_tag(tag: str):
    VALID_TAGS = [
        "18S_fungal_sequences",
        "Betacoronavirus",
        "ref_viroids_rep_genomes",
        "ref_viruses_rep_genomes",
        "refseq_select_rna",
        "refseq_select_prot",
        "refseq_protein",
        "refseq_rna",
        "16S_ribosomal_RNA",
        "ITS_RefSeq_Fungi",
        "28S_fungal_sequences",
        "ITS_eukaryote_sequences",
        "LSU_eukaryote_rRNA",
        "LSU_prokaryote_rRNA",
        "SSU_eukaryote_rRNA",
        "env_nt",
        "env_nr",
        "human_genome",
        "landmark",
        "mito",
        "mouse_genome",
        "nr",
        "nt_euk",
        "nt",
        "nt_others",
        "swissprot",
        "tsa_nr",
        "tsa_nt",
        "taxdb",
        "nt_prok",
        "nt_viruses",
        "pataa",
        "patnt",
        "pdbaa",
        "pdbnt",
        "ref_euk_rep_genomes",
        "ref_prok_rep_genomes",
    ]
    if tag not in VALID_TAGS:
        raise ValueError(f"{tag =} was inferred as Blast DB tag, which is not valid. {VALID_TAGS =}")


def get_blast_ref_tag():
    return parse_blast_tag(config["assembly__classification__blast"]["db_dir"])


if "blast" in config["assembly"]["classification"]:
    validate_blast_tag(get_blast_ref_tag())


### Global rule-set stuff #############################################################################################


BLAST_HEADER = "qseqid sacc staxid sscinames scomnames stitle pident evalue length mismatch gapopen qstart qend sstart send qlen slen"


def get_blast_binary():
    db_type = config["assembly__classification__blast"]["query_vs_db"]
    blast_binaries = {
        "nucleotide-nucleotide": "blastn",
        "nucleotide-protein": "blastx",
        "protein-nucleotide": "tblastn",
        "protein-protein": "blastp",
    }
    return blast_binaries[db_type]


def get_max_number_of_hits():
    return config["assembly__classification__blast"]["max_number_of_hits"]


def get_blast_db():
    blast_type = config["assembly__classification__blast"]["query_vs_db"].split("-")[1][0]
    reference_tag = get_blast_ref_tag()
    return multiext(
        os.path.join(config["assembly__classification__blast"]["db_dir"], f"{reference_tag}.{blast_type}"),
        "db",
        "ot",
        "tf",
        "to",
    )


def get_constraints():
    return {
        "reference": "|".join(get_reference_names()),
        "bam_step": "|".join(["original", "deduplication"]),
    }


def temp_mapping(output_file):
    if get_last_bam_step() == "original":
        return output_file
    return temp(output_file)


def get_last_bam_step():
    return "deduplication" if config["mapping"]["deduplication"] else "original"


def get_outputs_of_mapping():
    sample_names = get_sample_names()

    outputs = {}
    outputs["bams"] = expand(
        f"results/mapping/{{reference}}/{{sample}}.{get_last_bam_step()}.bam",
        sample=sample_names,
        reference=get_reference_names(),
    )

    if config["mapping"]["_generate_qualimap"]:
        outputs["qualimaps"] = expand(
            f"results/mapping/{{reference}}/bamqc/{{sample}}.{get_last_bam_step()}",
            sample=sample_names,
            reference=get_reference_names(),
        )

    return outputs


def get_assembly_fasta(sample: str):
    assembly_tool = config["assembly"]["assembly"]
    return f"results/assembly/{sample}/{assembly_tool}/contigs.fasta"


def get_outputs_of_assembly():
    if config["assembly"]["assembly"] == "":
        return {}

    outputs = {}
    sample_names = get_sample_names()

    assembly_outputs = []
    for sample in sample_names:
        assembly_outputs.append(get_assembly_fasta(sample))
    outputs["assembly"] = assembly_outputs

    if "quast" in config["assembly"]["report"]:
        outputs["quast"] = expand(
            "results/assembly/{sample}/spades/QUAST/report.pdf",
            sample=sample_names,
        )
    if "bandage" in config["assembly"]["report"]:
        outputs["bandage"] = expand(
            "results/assembly/{sample}/spades/bandage/bandage.{ext}",
            sample=sample_names,
            ext=["info", "svg"],
        )

    if "blast" in config["assembly"]["classification"]:
        outputs["blast"] = expand(
            "results/classification/{sample}/spades/blast/summary.html",
            sample=sample_names,
        )

    return outputs


def get_outputs():
    return get_outputs_of_mapping() | get_reads_outputs() | get_outputs_of_assembly()


def get_standalone_outputs():
    # outputs that will be produced if the module is run as a standalone workflow, not as a part of a larger workflow
    return {
        "multiqc_report": expand("results/_aggregation/multiqc_{reference}.html", reference=get_reference_names()),
    }


def infer_bwa_index_for_mapping(wildcards):
    return multiext(
        os.path.join(get_reference_dir_for_name(wildcards.reference), "bwa_index", wildcards.reference),
        ".amb",
        ".ann",
        ".bwt",
        ".pac",
        ".sa",
    )


def infer_hisat2_index_for_mapping(wildcards):
    return multiext(
        os.path.join(get_reference_dir_for_name(wildcards.reference), "hisat2_index", wildcards.reference),
        ".1.ht2l",
        ".2.ht2l",
        ".3.ht2l",
        ".4.ht2l",
        ".5.ht2l",
        ".6.ht2l",
        ".7.ht2l",
        ".8.ht2l",
    )


def infer_dragmap_index_for_mapping(wildcards):
    return multiext(
        os.path.join(get_reference_dir_for_name(wildcards.reference), "dragmap_index", wildcards.reference),
        "hash_table.cfg",
        "hash_table.cfg.bin",
        "hash_table.cmp",
        "hash_table_stats.txt",
        "reference.bin",
        "ref_index.bin",
        "repeat_mask.bin",
        "str_table.bin",
    )


def infer_minimap2_index_for_mapping(wildcards):
    return os.path.join(get_reference_dir_for_name(wildcards.reference), f"{wildcards.reference}.mmi")


def infer_final_bam(wildcards):
    return get_input_bam_for_sample_and_ref(wildcards.sample, wildcards.reference)


def infer_final_bai(wildcards):
    return get_input_bai_for_sample_and_ref(wildcards.sample, wildcards.reference)


def infer_reference_fasta(wildcards):
    return reference_dict[wildcards.reference]


### Contract for other workflows ######################################################################################


def get_input_bam_for_sample_and_ref(sample: str, reference: str):
    return f"results/mapping/{reference}/{sample}.{get_last_bam_step()}.bam"


def get_input_bai_for_sample_and_ref(sample: str, reference: str):
    return f"results/mapping/{reference}/{sample}.{get_last_bam_step()}.bam.bai"


def get_multiqc_inputs(reference: str):
    outs = get_multiqc_inputs_for_reads()

    if config["mapping"]["deduplication"] in ["picard", "samtools"]:
        outs["picard_dedup"] = expand(
            f"results/mapping/{reference}/{{sample}}.deduplication.stats",
            sample=get_sample_names(),
        )

    if config["mapping"]["_generate_qualimap"]:
        outs["qualimaps"] = expand(
            f"results/mapping/{reference}/bamqc/{{sample}}.{get_last_bam_step()}",
            sample=get_sample_names(),
        )

    outs["samtools_stats"] = expand(
        f"results/mapping/{reference}/{{sample}}.{{step}}.samtools_stats",
        step=get_last_bam_step(),
        sample=get_sample_names(),
    )
    return outs


def infer_multiqc_inputs_for_reference(wildcards):
    return get_multiqc_inputs(wildcards.reference)


### Parameter parsing from config #####################################################################################


def get_spades_params():
    mode = (
        ""
        if config["assembly__assembly__spades"]["mode"] == "standard"
        else f'--{config["assembly__assembly__spades"]["mode"]}'
    )
    careful = "--careful" if config["assembly__assembly__spades"]["careful"] else ""
    if mode and careful:
        return f"{mode} {careful}"
    return mode + careful


def get_quast_params():
    mincontig_param = "--min-contig {val}".format(val=config["assembly__report__quast"]["min_contig_length"])
    if config["assembly__report__quast"]["extra"]:
        return f'{mincontig_param} {config["assembly__report__quast"]["extra"]}'
    return mincontig_param


### Resource handling #################################################################################################


def get_mem_mb_for_assembly(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["assembly__assembly_mem_mb"] * attempt)


def get_mem_mb_for_classification(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["assembly__classification_mem_mb"] * attempt)


def get_mem_mb_for_deduplication(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["mapping__deduplication_mem_mb"] * attempt)


def get_mem_mb_for_qualimap(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["mapping__qualimap_mem_mb"] * attempt)


def get_mem_mb_for_mapping(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["mapping__mapping_mem_mb"] * attempt)


def get_mem_mb_for_indexing(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["mapping__indexing_mem_mb"] * attempt)


def get_threads_for_assembly():
    return min(config["threads"]["assembly__assembly"], config["max_threads"])


def get_threads_for_classification():
    return min(config["threads"]["assembly__classification"], config["max_threads"])


def get_threads_for_report():
    return min(config["threads"]["assembly__report"], config["max_threads"])
