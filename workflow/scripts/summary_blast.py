import csv
import sys

import pandas as pd


def summary_blast(
    blast_tsv: str,
    min_query_coverage: float,
    max_target_seqs: int,
    output_tsv: str,
):
    seqids, info_list = [], []

    blast = pd.read_csv(blast_tsv, sep="\t", index_col=None, encoding="utf-8", quoting=csv.QUOTE_NONE)

    for seqid, aligns in blast.groupby("qseqid"):
        aligns["qcov"] = ((aligns["qend"] - aligns["qstart"]) / aligns["qlen"]).abs()
        aligns = aligns[aligns["qcov"] >= min_query_coverage]
        aligns = aligns.iloc[:max_target_seqs]

        if len(aligns) == 0:
            continue

        def create_link(align):
            return '%3.2f%% - <a href="https://www.ncbi.nlm.nih.gov/nuccore/%s", title="%s">%s</a>' % (
                align["qcov"] * 100,
                align["sacc"],
                align["sscinames"],
                align["stitle"],
            )

        seqids.append(seqid)
        info_list.append(
            {
                "Homologue title": ";".join(aligns["stitle"]),
                "Homologue accession": ";".join(aligns["stitle"]),
                "Homologue link": "<br />".join(aligns.apply(create_link, axis=1)),
            }
        )

    infos = pd.DataFrame(info_list, index=seqids)
    infos.to_csv(output_tsv, sep="\t")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")

    summary_blast(
        snakemake.input.tsv,
        snakemake.params.min_query_coverage,
        snakemake.params.max_target_seqs,
        snakemake.output.tsv,
    )
