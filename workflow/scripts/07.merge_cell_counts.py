#!/usr/bin/env python3
import glob
import sys
from pathlib import Path

import pandas as pd


def main():
    if len(sys.argv) not in (2, 4):
        print(
            "Usage: 07.merge_cell_counts.py <dataset> [input_dir output_file]",
            file=sys.stderr,
        )
        sys.exit(1)

    dataset = sys.argv[1]
    input_dir = Path(sys.argv[2]) if len(sys.argv) == 4 else Path(f"results/06.Count_cells/{dataset}")
    output_file = (
        Path(sys.argv[3])
        if len(sys.argv) == 4
        else Path("results/07.merge_cell_counts") / f"{dataset}_all_samples_cell_counts.tsv"
    )
    output_file.parent.mkdir(parents=True, exist_ok=True)

    all_files = sorted(glob.glob(str(input_dir / f"{dataset}_*_cell_counts.tsv")))
    if not all_files:
        raise FileNotFoundError(f"No cell count files found in {input_dir} for dataset {dataset}")

    all_counts = [pd.read_csv(f, sep="\t") for f in all_files]
    merged = pd.concat(all_counts, axis=0, ignore_index=True)

    # Add a total row for numeric columns.
    sum_row = merged.sum(numeric_only=True)
    sum_row.name = "Total"
    merged = pd.concat([merged, pd.DataFrame(sum_row).T], axis=0, ignore_index=False)

    merged.to_csv(output_file, sep="\t", index=True)
    print(f"Merged {len(all_files)} files into {output_file}")


if __name__ == "__main__":
    main()
