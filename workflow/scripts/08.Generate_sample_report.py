#!/usr/bin/env python3
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import image as mpimg
from matplotlib.backends.backend_pdf import PdfPages


def main():
    if len(sys.argv) != 8:
        print(
            "Usage: 08.Generate_sample_report.py <dataset> <sample> <out_pdf> "
            "<knee_plot.png> <scrublet_hist.png> <qc_distributions.png> <cell_counts.tsv>",
            file=sys.stderr,
        )
        sys.exit(1)

    dataset = sys.argv[1]
    sample = sys.argv[2]
    out_pdf = sys.argv[3]
    knee_plot_file = sys.argv[4]
    scrublet_hist_file = sys.argv[5]
    qc_dist_file = sys.argv[6]
    cell_counts_file = sys.argv[7]

    Path(out_pdf).parent.mkdir(parents=True, exist_ok=True)

    for path in [knee_plot_file, scrublet_hist_file, qc_dist_file, cell_counts_file]:
        if not Path(path).exists():
            raise FileNotFoundError(path)

    df = pd.read_csv(cell_counts_file, sep="\t")

    with PdfPages(out_pdf) as pdf:
        fig, axes = plt.subplots(
            4,
            1,
            figsize=(8, 10),
            gridspec_kw={"height_ratios": [3, 1.5, 3, 0.5]},
            dpi=600,
        )

        img = mpimg.imread(knee_plot_file)
        axes[0].imshow(img)
        axes[0].axis("off")
        axes[0].set_title(f"{dataset} {sample} - DropletUtils Knee Plot", fontsize=14)

        img = mpimg.imread(scrublet_hist_file)
        axes[1].imshow(img)
        axes[1].axis("off")
        axes[1].set_title(f"{dataset} {sample} - Scrublet Histogram", fontsize=14)

        img = mpimg.imread(qc_dist_file)
        axes[2].imshow(img)
        axes[2].axis("off")
        axes[2].set_title(f"{dataset} {sample} - scAutoQC QC Distributions", fontsize=14)

        axes[3].axis("off")
        axes[3].axis("tight")
        table = axes[3].table(
            cellText=df.values,
            colLabels=df.columns,
            loc="center",
            cellLoc="center",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 1.5)
        axes[3].set_title(f"{dataset} {sample} - Cell Counts", fontsize=14)

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()


if __name__ == "__main__":
    main()
