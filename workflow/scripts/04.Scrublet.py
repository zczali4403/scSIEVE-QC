#!/usr/bin/env python3
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scrublet as scr


def main():
    if len(sys.argv) not in (4, 6):
        print(
            "Usage: 04.Scrublet.py <input.h5ad> <output_filtered.h5ad> "
            "<output_hist.png> [expected_doublet_rate] [threshold]",
            file=sys.stderr,
        )
        sys.exit(1)

    in_file = sys.argv[1]
    out_file = sys.argv[2]
    out_png = sys.argv[3]
    expected_doublet_rate = float(sys.argv[4]) if len(sys.argv) >= 5 else 0.06
    threshold = float(sys.argv[5]) if len(sys.argv) >= 6 else 0.4

    Path(out_file).parent.mkdir(parents=True, exist_ok=True)
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(in_file)

    scrub = scr.Scrublet(adata.X, expected_doublet_rate=expected_doublet_rate)
    doublet_scores, _ = scrub.scrub_doublets(
        min_counts=2,
        min_cells=3,
        min_gene_variability_pctl=85,
        n_prin_comps=30,
    )

    # Save histogram. If Scrublet did not infer a threshold, draw a fallback plot.
    if hasattr(scrub, "threshold_") and scrub.threshold_ is not None:
        scrub.plot_histogram()
    else:
        bins = np.linspace(0, 1, 50)
        fig, axs = plt.subplots(1, 2, figsize=(8, 3))
        ax = axs[0]
        ax.hist(doublet_scores, bins=bins, color="gray", linewidth=0, density=True)
        ax.set_yscale("log")
        ax.set_title("Observed transcriptomes")
        ax.set_xlabel("Doublet score")
        ax.set_ylabel("Prob. density")

        ax = axs[1]
        ax.hist(scrub.doublet_scores_sim_, bins=bins, color="gray", linewidth=0, density=True)
        ax.set_yscale("linear")
        ax.set_title("Simulated doublets")
        ax.set_xlabel("Doublet score")
        ax.set_ylabel("Prob. density")
        fig.tight_layout()

    plt.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.close()

    scrub.call_doublets(threshold=threshold)
    adata.obs["doublet_score"] = doublet_scores
    adata.obs["doublets"] = doublet_scores > threshold

    adata = adata[~adata.obs["doublets"]].copy()
    adata.write(out_file)


if __name__ == "__main__":
    main()
