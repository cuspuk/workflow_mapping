from snakemake.utils import validate


configfile: "config/config.yaml"


validate(config, "../schemas/config.schema.yaml", set_default=False)


### Layer for adapting other workflows  ###############################################################################


def get_fastq_for_mapping(wildcards):
    return reads_workflow.get_final_fastq_for_sample(wildcards.sample)


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


# def sample_has_enough_mapped_reads(sample_name: str):
#     with checkpoints.checkpoint__min_mapped_reads_count.get(sample=sample_name).output[0].open() as f:
#         num = int(f.read().strip())
#     return num >= 1

### Global rule-set stuff #############################################################################################


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


def get_outputs():
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

    outputs = outputs | get_reads_outputs()
    return outputs


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
    return os.path.join(get_reference_dir_for_name(wildcards.reference), f"{wildcards.fasta}.mmi")


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


### Resource handling #################################################################################################


def get_mem_mb_for_deduplication(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["mapping__deduplication_mem_mb"] * attempt)


def get_mem_mb_for_qualimap(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["mapping__qualimap_mem_mb"] * attempt)


def get_mem_mb_for_mapping(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["mapping__mapping_mem_mb"] * attempt)


def get_mem_mb_for_indexing(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["mapping__indexing_mem_mb"] * attempt)
