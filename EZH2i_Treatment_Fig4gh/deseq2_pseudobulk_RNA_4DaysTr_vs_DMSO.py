#!/usr/bin/env python3
"""Reproduce the All_cells / selected / 4 Days Tr. versus DMSO DESeq2 outputs."""
import argparse
from pathlib import Path
import subprocess
import tempfile

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text

PREFIX = "01__4_Days_Tr_vs_DMSO__DESeq2"
R_CODE = r'''
args <- commandArgs(TRUE)
suppressPackageStartupMessages(library(DESeq2))
x <- read.delim(args[1], check.names=FALSE, row.names=1)
counts <- as.matrix(x[, 1:6])
storage.mode(counts) <- "integer"
group <- factor(rep(c("ref", "alt"), each=3), levels=c("ref", "alt"))
dds <- DESeqDataSetFromMatrix(counts,
    data.frame(group=group, row.names=colnames(counts)), design=~group)
dds <- DESeq(dds)
res <- results(dds, contrast=c("group", "alt", "ref"))
out <- data.frame(gene=rownames(res), logFC=res$log2FoldChange,
    logCPM=log2(res$baseMean + 0.5), PValue=res$pvalue, FDR=res$padj)
write.table(out, args[2], sep="\t", row.names=FALSE, quote=FALSE, na="NA")
cat("R:", as.character(getRversion()), "; DESeq2:",
    as.character(packageVersion("DESeq2")), "\n")
'''


def pseudobulk_from_h5ad(path):
    a = ad.read_h5ad(path)
    conditions = a.obs["Condition"].astype(str).to_numpy()
    rng = np.random.default_rng(12345)
    # Replay only metadata permutations from the four earlier contrasts to
    # preserve the original selected/01 grouping; no other expression is used.
    for ref, alt in [("DMSO", "4 Days Tr."), ("4 Days Tr.", "2 Days Rec."),
                     ("2 Days Rec.", "4 Days Rec."), ("4 Days Rec.", "8 Days Rec.")]:
        labels = conditions[np.isin(conditions, [ref, alt])]
        for condition in (ref, alt):
            rng.permutation(np.flatnonzero(labels == condition))

    names = a.var_names.astype(str).str.upper()
    keep = ~(names.str.contains("MALAT1") |
             names.str.startswith(("MT-", "MT_", "MT.")))
    a = a[np.isin(conditions, ["DMSO", "4 Days Tr."]), keep].copy()
    if not a.var_names.is_unique:
        raise ValueError("Gene names must be unique")
    X = a.layers["counts"]
    labels = a.obs["Condition"].astype(str).to_numpy()
    units, columns, cell_counts = [], [], []
    for label, condition, slug in [("ref", "DMSO", "DMSO"),
                                   ("alt", "4 Days Tr.", "4_Days_Tr")]:
        idx = np.flatnonzero(labels == condition)
        if len(idx) < 25:
            raise ValueError(f"{condition} needs at least 25 cells; found {len(idx)}")
        cell_counts.append(len(idx))
        for rep, chunk in enumerate(np.array_split(rng.permutation(idx), 3), 1):
            units.append(np.rint(np.asarray(X[chunk].sum(axis=0)).ravel()).astype(np.int64))
            columns.append(f"{label}__{slug}__rep{rep}")
    table = pd.DataFrame(np.column_stack(units), index=a.var_names, columns=columns)
    table.index.name = "gene"
    table["n_cells_detected"] = np.asarray((X > 0).sum(axis=0)).ravel().astype(float)
    print(f"DMSO: {cell_counts[0]} cells; 4 Days Tr.: {cell_counts[1]} cells; "
          f"{a.n_vars} genes; using layers['counts']", flush=True)
    return table, cell_counts


def volcano(df, outbase):
    d = df.dropna(subset=["gene", "logFC", "FDR"]).copy()
    d["mlog10FDR"] = -np.log10(d.FDR.clip(lower=np.finfo(float).tiny))
    up = (d.FDR < 0.01) & (d.logFC >= 1)
    down = (d.FDR < 0.01) & (d.logFC <= -1)
    sizes = 6 + 54 * np.sqrt(np.clip(d.pct_cells_detected.to_numpy(), 0, 1))
    fig, ax = plt.subplots(figsize=(5.8, 5.8), dpi=220)
    for mask, color, alpha, label in [
        (~(up | down), "#c7c7c7", 0.45, "NS"),
        (down, "#3b82f6", 0.85, "Higher in DMSO"),
        (up, "#ef4444", 0.85, "Higher in 4 Days Tr."),
    ]:
        ax.scatter(d.loc[mask, "logFC"], d.loc[mask, "mlog10FDR"],
                   s=sizes[mask], c=color, alpha=alpha, linewidths=0, label=label)
    for x in (1, -1):
        ax.axvline(x, color="black", ls="--", lw=1, alpha=0.5)
    ax.axhline(2, color="black", ls="--", lw=1, alpha=0.5)
    ax.set_xlabel("log2 fold-change (4 Days Tr. − DMSO)")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title(
        "All_cells | selected | DESeq2\nlogFC = 4 Days Tr. − DMSO\n"
        f"Up in ALT: {int(up.sum())}   Up in REF: {int(down.sum())}   "
        "FDR < 0.01, |log2FC| ≥ 1.0", fontsize=9)
    ax.legend(frameon=False, fontsize=7, loc="best")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.grid(True, axis="y", ls=":", lw=0.6, alpha=0.6)
    top = (d.loc[up | down].assign(abs_logFC=lambda x: x.logFC.abs())
           .sort_values(["FDR", "abs_logFC"], ascending=[True, False]).head(80))
    texts = [ax.text(r.logFC, r.mlog10FDR, r.gene, fontsize=6,
                     ha="center", va="bottom") for r in top.itertuples()]
    if texts:
        adjust_text(texts, ax=ax, expand_points=(1.2, 1.4),
                    expand_text=(1.1, 1.2),
                    arrowprops=dict(arrowstyle="-", lw=0.5, color="black", alpha=0.6))
    fig.tight_layout()
    fig.savefig(str(outbase) + "__volcano.pdf")
    fig.savefig(str(outbase) + "__volcano.png", dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("rna.h5ad"),
                        help="RNA AnnData file (default: ./rna.h5ad)")
    parser.add_argument("--outdir", type=Path,
                        default=Path("results/DESeq2_RNA_4DaysTr_vs_DMSO"),
                        help="Output folder (default: ./results/DESeq2_RNA_4DaysTr_vs_DMSO)")
    parser.add_argument("--rscript", default="Rscript")
    args = parser.parse_args()
    x, (n_ref, n_alt) = pseudobulk_from_h5ad(args.input)
    with tempfile.TemporaryDirectory() as tmp:
        rfile, result = Path(tmp) / "run.R", Path(tmp) / "results.tsv"
        counts_file = Path(tmp) / "counts.tsv"
        x.to_csv(counts_file, sep="\t")
        rfile.write_text(R_CODE)
        subprocess.run([args.rscript, str(rfile), str(counts_file),
                        str(result)], check=True)
        df = pd.read_csv(result, sep="\t", keep_default_na=False, na_values=["NA"])
    assert df.gene.tolist() == x.index.tolist(), "Input/output gene order differs"
    df["n_cells_detected"] = x.n_cells_detected.to_numpy()
    df["pct_cells_detected"] = df.n_cells_detected / (n_ref + n_alt)
    metadata = dict(method="DESeq2", epistate_stratum="All_cells",
                    epistate_values="ALL", contrast_set="selected", contrast_index=1,
                    ref_condition="DMSO", alt_condition="4 Days Tr.",
                    logFC_direction="4 Days Tr. - DMSO", n_cells_ref=n_ref,
                    n_cells_alt=n_alt, n_pseudobulk_ref=3, n_pseudobulk_alt=3)
    for i, (key, value) in enumerate(metadata.items()):
        df.insert(i, key, value)
    args.outdir.mkdir(parents=True, exist_ok=True)
    outbase = args.outdir / PREFIX
    df.to_csv(str(outbase) + ".csv", index=False)
    np.random.seed(12345)
    volcano(df, outbase)
    print(f"Wrote CSV, volcano PDF and volcano PNG to {args.outdir}")


if __name__ == "__main__":
    main()
