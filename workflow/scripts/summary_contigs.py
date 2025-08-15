import os
import shutil
import sys
import zlib

import numpy as np
import pandas as pd
from Bio import SeqIO, SeqUtils


def split_contigs(input_fasta: str, params_fasta: str, output_fasta_dir: str, output_seqinfo: str):
    os.makedirs(output_fasta_dir, exist_ok=True)
    shutil.copyfile(input_fasta, params_fasta)

    seqids, info_list = [], []
    for seq in SeqIO.parse(input_fasta, "fasta"):
        print(f"Processing {seq.id}", file=sys.stderr)
        seq_path = os.path.join(output_fasta_dir, seq.id)
        with open(seq_path, "w") as out:
            SeqIO.write(seq, out, "fasta")

        seqids.append(seq.id)
        contig_info = {
            "Length": len(seq),
            "GC content": str(np.round(SeqUtils.GC(seq.seq), 1)),
            "Compress ratio": len(zlib.compress(bytes(str(seq.seq), "utf-8"))) / len(seq.seq),
        }
        info_list.append(contig_info)

    infos = pd.DataFrame(info_list, index=seqids)
    infos.to_csv(output_seqinfo, sep="\t")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")

    split_contigs(
        snakemake.input.contigs,
        snakemake.params.fasta,
        snakemake.output.fasta_dir,
        snakemake.output.seqinfo,
    )
