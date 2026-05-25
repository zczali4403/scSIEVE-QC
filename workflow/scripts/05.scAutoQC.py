#!/usr/bin/env python3
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import sctk
import seaborn as sns


def main():
    if len(sys.argv) not in (5, 7):
        print(
            "Usage: 05.scAutoQC.py <input.h5ad> <out_dir> <dataset> <sample> "
            "[min_genes] [n_counts_max_quantile]",
            file=sys.stderr,
        )
        sys.exit(1)

    in_file = sys.argv[1]
    out_dir = sys.argv[2]
    dataset = sys.argv[3]
    sample = sys.argv[4]
    min_genes = int(sys.argv[5]) if len(sys.argv) >= 6 else 200
    n_counts_max_quantile = float(sys.argv[6]) if len(sys.argv) >= 7 else 0.95

    Path(out_dir).mkdir(parents=True, exist_ok=True)
    out_file_prefix = f"{out_dir}/{dataset}_{sample}"

    sc.settings.set_figure_params(dpi=80)

    adata = sc.read_h5ad(in_file)
    print(adata)

    # Compute QC metrics.
    sctk.calculate_qc(adata)
    print(adata)

    # Configure cell-wise QC thresholds.
    metrics = sctk.default_metric_params_df.loc[
        ["n_counts", "n_genes", "percent_mito", "percent_ribo"], :
    ].copy()
    metrics.loc["n_genes", "min"] = min_genes
    metrics.loc["n_counts", "min"] = 0
    metrics.loc["n_counts", "max"] = adata.obs["n_counts"].quantile(n_counts_max_quantile)
    metrics.loc["n_genes", "side"] = "both"
    metrics.loc["n_counts", "side"] = "both"

    sctk.cellwise_qc(adata, metrics=metrics)
    print(adata)
    print("cell_passed_qc:", adata.obs["cell_passed_qc"].sum())

    # QC clustering.
    metrics_list = ["log1p_n_counts", "log1p_n_genes", "percent_mito", "percent_ribo"]
    sctk.generate_qc_clusters(adata, metrics=metrics_list)
    print(adata)

    sc.pl.embedding(
        adata,
        "X_umap_qc",
        color=["qc_cluster", "log1p_n_counts"],
        color_map="OrRd",
        ncols=2,
        wspace=0.4,
        show=False,
    )
    plt.savefig(f"{out_file_prefix}_umap_qc_cluster.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Cluster-wise QC.
    sctk.clusterwise_qc(adata)
    print(adata)
    print(adata.obs.cluster_passed_qc.value_counts())

    for col in ["cell_passed_qc", "cluster_passed_qc"]:
        adata.obs[col + "_int"] = adata.obs[col].astype(int)

    sc.pl.embedding(
        adata,
        "X_umap_qc",
        color=["cell_passed_qc_int", "cluster_passed_qc_int"],
        show=False,
    )
    plt.savefig(f"{out_file_prefix}_umap_qc_cluster_passQC.png", dpi=300, bbox_inches="tight")
    plt.close()

    # QC distribution plots.
    all_counts = adata.obs["n_counts"]
    all_genes = adata.obs["n_genes"]
    qc_pass = adata.obs["cluster_passed_qc_int"] == 1
    qc_counts = adata.obs.loc[qc_pass, "n_counts"]
    qc_genes = adata.obs.loc[qc_pass, "n_genes"]

    ranges = adata.uns["scautoqc_ranges"]
    counts_range = ranges.loc["n_counts"].tolist()
    genes_range = ranges.loc["n_genes"].tolist()
    print("n_counts range:", counts_range)
    print("n_genes range:", genes_range)

    x_counts_min = min(all_counts.min(), qc_counts.min())
    x_counts_max = max(all_counts.max(), qc_counts.max())
    x_genes_min = min(all_genes.min(), qc_genes.min())
    x_genes_max = max(all_genes.max(), qc_genes.max())

    bins_counts = np.linspace(x_counts_min, x_counts_max, 50)
    bins_genes = np.linspace(x_genes_min, x_genes_max, 50)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    sns.histplot(all_counts, bins=bins_counts, kde=False, ax=axes[0, 0], color="gray")
    axes[0, 0].axvspan(counts_range[0], counts_range[1], color="red", alpha=0.2)
    axes[0, 0].set_title("All cells - n_counts")

    sns.histplot(all_genes, bins=bins_genes, kde=False, ax=axes[0, 1], color="gray")
    axes[0, 1].axvspan(genes_range[0], genes_range[1], color="red", alpha=0.2)
    axes[0, 1].set_title("All cells - n_genes")

    sns.histplot(qc_counts, bins=bins_counts, kde=False, ax=axes[1, 0], color="green")
    axes[1, 0].set_title("QC passed - n_counts")

    sns.histplot(qc_genes, bins=bins_genes, kde=False, ax=axes[1, 1], color="green")
    axes[1, 1].set_title("QC passed - n_genes")

    ylim_counts = max(axes[0, 0].get_ylim()[1], axes[1, 0].get_ylim()[1])
    ylim_genes = max(axes[0, 1].get_ylim()[1], axes[1, 1].get_ylim()[1])
    axes[0, 0].set_ylim(0, ylim_counts)
    axes[1, 0].set_ylim(0, ylim_counts)
    axes[0, 1].set_ylim(0, ylim_genes)
    axes[1, 1].set_ylim(0, ylim_genes)

    plt.tight_layout()
    plt.savefig(f"{out_file_prefix}_qc_distributions.png", dpi=600)
    plt.close()

    # Save filtered AnnData.
    qc_pass_cells = adata.obs["cluster_passed_qc"]
    adata_filtered = adata[qc_pass_cells].copy()
    adata_filtered.write(f"{out_file_prefix}_scAutoQC_adata.h5ad")


if __name__ == "__main__":
    main()
