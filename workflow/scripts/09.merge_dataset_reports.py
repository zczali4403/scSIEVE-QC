#!/usr/bin/env python3
import glob
import os
import sys
from pathlib import Path

from PyPDF2 import PdfReader, PdfWriter


def main():
    if len(sys.argv) < 3:
        print(
            "Usage: 09.merge_dataset_reports.py <dataset_name> <output_pdf> [input_pdf ...]",
            file=sys.stderr,
        )
        sys.exit(1)

    dataset = sys.argv[1]
    output_file = Path(sys.argv[2])
    input_files = sys.argv[3:]

    # Backward-compatible fallback if explicit input files were not supplied.
    if not input_files:
        input_files = sorted(glob.glob(f"results/08.Generate_sample_report/{dataset}/*_report.pdf"))

    input_files = sorted(input_files)
    if not input_files:
        raise FileNotFoundError(f"No input PDF files found for dataset {dataset}")

    output_file.parent.mkdir(parents=True, exist_ok=True)
    writer = PdfWriter()

    for pdf_file in input_files:
        reader = PdfReader(pdf_file)
        for page in reader.pages:
            writer.add_page(page)

    with open(output_file, "wb") as f:
        writer.write(f)

    print(f"Merged {len(input_files)} PDFs into {output_file}")


if __name__ == "__main__":
    main()
