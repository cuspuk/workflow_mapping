import os
import sys

import numpy as np
import pandas as pd

pd.set_option("display.max_colwidth", 1000000)
pd.options.display.float_format = "{:.3f}".format


def load_attrs(attr_file: str):
    table = pd.read_csv(attr_file, sep="\t", index_col=0)
    return table


def build_report(
    *,
    html_template: str,
    fasta: str,
    blast_file: str,
    seqinfo: str,
    max_query_seqs: int,
    output_tsv: str,
    output_html: str,
    fasta_dir: str,
    html_columns_req: list[str],
    sort_by: str,
    sort_how: str,
    seqs_per_page: int,
):
    attr_files = [blast_file] + [seqinfo]
    attr_names = [os.path.basename(f)[:-4] for f in attr_files]
    attr_tables = [load_attrs(attr_file) for attr_file in attr_files]

    attrs, last_name = attr_tables[0], attr_names[0]
    for attr_name, attr_table in zip(attr_names[1:], attr_tables[1:]):
        attrs = attrs.merge(attr_table, how="outer", left_index=True, right_index=True)
        last_name = attr_name

    attrs.sort_values(by="Length", ascending=False, inplace=True)
    if max_query_seqs:
        attrs = attrs.iloc[:max_query_seqs]

    tsv_columns = [attr for attr in attrs.columns if " link" not in attr]
    attrs[tsv_columns].to_csv(output_tsv, sep="\t")

    TEMPLATE = open(html_template).read()
    with open(output_html, "w") as out:
        attrs["Sequence"] = [
            '<a href="%s/%s.fa">%s</a>' % (fasta_dir, seqid, i + 1) for i, seqid in enumerate(attrs.index)
        ]
        attrs.index = np.arange(1, len(attrs) + 1)

        html_columns = attrs.columns
        if html_columns_req:
            html_columns = []
            for col_prefix in html_columns_req:
                html_columns.extend(sorted([col for col in attrs.columns if col.startswith(col_prefix)]))

        html_attrs = attrs[html_columns]
        html_table = html_attrs.fillna("").to_html(escape=False, index=False).replace("<table ", '<table id="data" ')

        sort_index = list(html_attrs.columns).index(sort_by)
        out.write(
            TEMPLATE.format(
                attrs=html_table, seq_fasta=fasta, sort_index=sort_index, sort_how=sort_how, seqs_per_page=seqs_per_page
            )
        )


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")

    build_report(
        html_template=snakemake.input.template,
        fasta=snakemake.params.fasta,
        blast_file=snakemake.input.attrs_blast,
        seqinfo=snakemake.input.attrs_seqinfo,
        max_query_seqs=snakemake.params.max_query_seqs,
        seqs_per_page=snakemake.params.seqs_per_page,
        sort_by=snakemake.params.sort_by,
        sort_how=snakemake.params.sort_how,
        output_html=snakemake.output.html,
        output_tsv=snakemake.output.tsv,
        fasta_dir=snakemake.params.fasta_dir,
        html_columns_req=snakemake.params.html_columns,
    )
