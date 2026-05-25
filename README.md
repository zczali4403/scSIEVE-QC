# scSIEVE-QC

**scSIEVE-QC** stands for **single-cell Snakemake-Integrated Empty-droplet filtering, doublet Verification, and Expression-quality Evaluation for Quality Control**.

scSIEVE-QC is a Snakemake workflow for scalable single-cell RNA-seq quality control. It integrates:

- **DropletUtils** for empty-droplet filtering
- **Scrublet** for doublet detection and removal
- **scAutoQC / sctk** for automated cell- and cluster-level quality evaluation
- **Snakemake** for reproducible workflow orchestration

The workflow takes per-sample raw/unfiltered `h5ad` files as input and produces filtered objects, QC plots, cell-count summaries, per-sample PDF reports, and dataset-level merged reports.

---

## Workflow overview

For each sample, scSIEVE-QC runs:

1. **h5ad to Seurat**  
   Convert the raw input `h5ad` file to a Seurat RDS object.

2. **Seurat to SingleCellExperiment**  
   Convert Seurat to `SingleCellExperiment` for DropletUtils.

3. **DropletUtils empty-droplet filtering**  
   Estimate barcode-rank statistics, run `emptyDrops`, generate a knee plot, and retain cells using the workflow's UMI thresholding strategy.

4. **Seurat to h5ad**  
   Convert the filtered Seurat object back to `h5ad`.

5. **Scrublet doublet removal**  
   Compute doublet scores, save a Scrublet histogram, and remove predicted doublets.

6. **scAutoQC expression-quality evaluation**  
   Compute QC metrics, run cell-wise and cluster-wise QC, generate QC UMAPs and distribution plots, and save the final filtered `h5ad`.

7. **Cell-count summary**  
   Count retained cells after DropletUtils, Scrublet, and scAutoQC.

8. **Sample report**  
   Generate a per-sample PDF report.

9. **Dataset report**  
   Merge all sample-level reports within each dataset.

---

## Repository structure

```text
scSIEVE-QC/
├── README.md
├── README_zh.md
├── LICENSE
├── .gitignore
├── config/
│   ├── config.yaml
│   └── sample.csv
├── examples/
│   └── sample.csv
└── workflow/
    ├── Snakefile
    ├── envs/
    │   ├── sceasy.yaml
    │   └── scautoqc.yaml
    └── scripts/
        ├── 00.from_h5ad_to_seurat.r
        ├── 01.from_seurat_to_sce.r
        ├── 02.DropletUtils.r
        ├── 03.from_seurat_to_h5ad.r
        ├── 04.Scrublet.py
        ├── 05.scAutoQC.py
        ├── 06.Count_cells.py
        ├── 07.merge_cell_counts.py
        ├── 08.Generate_sample_report.py
        └── 09.merge_dataset_reports.py
```

Large input data and generated results are intentionally excluded from the repository.

---

## Input

Prepare a sample table with three required columns:

```csv
dataset,sample,file
Example_2024,sample1,/path/to/sample1_counts_unfiltered.h5ad
Example_2024,sample2,/path/to/sample2_counts_unfiltered.h5ad
```

Column meaning:

| Column | Description |
|---|---|
| `dataset` | Dataset name. Samples with the same dataset name are merged in dataset-level summaries. |
| `sample` | Sample identifier. |
| `file` | Path to the raw/unfiltered input `h5ad` file. Absolute paths are recommended. |

The default sample table is:

```text
config/sample.csv
```

You can edit `config/config.yaml` to point to another sample table.

---

## Configuration

Main configuration file:

```text
config/config.yaml
```

Important options:

```yaml
sample_table: "config/sample.csv"
outdir: "results"

envs:
  sceasy: "workflow/envs/sceasy.yaml"
  scautoqc: "workflow/envs/scautoqc.yaml"

scrublet:
  expected_doublet_rate: 0.06
  threshold: 0.4

scautoqc:
  min_genes: 200
  n_counts_max_quantile: 0.95
```

---

## Installation

Install Snakemake first, for example:

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake
```

Then run the workflow with Snakemake-managed conda environments:

```bash
snakemake -s workflow/Snakefile --use-conda --cores 8
```

### Notes on `sceasy`

The R package `sceasy` may not be available from all conda channels. If the environment created from `workflow/envs/sceasy.yaml` does not contain `sceasy`, install it in that environment following the upstream instructions, for example:

```bash
Rscript -e 'remotes::install_github("cellgeni/sceasy")'
```

---

## Usage

### 1. Edit the sample table

```bash
cp examples/sample.csv config/sample.csv
# edit config/sample.csv with your dataset/sample/file paths
```

### 2. Dry run

```bash
snakemake -s workflow/Snakefile \
  --configfile config/config.yaml \
  --use-conda \
  --cores 8 \
  -n
```

### 3. Run

```bash
snakemake -s workflow/Snakefile \
  --configfile config/config.yaml \
  --use-conda \
  --cores 8 \
  --rerun-incomplete
```

### 4. Print shell commands while running

```bash
snakemake -s workflow/Snakefile \
  --configfile config/config.yaml \
  --use-conda \
  --cores 8 \
  --rerun-incomplete \
  -p
```

### 5. Unlock after an interrupted run

```bash
snakemake -s workflow/Snakefile --unlock
```

---

## Output

All outputs are written under the configured `outdir`, default:

```text
results/
```

Main outputs:

```text
results/00.from_h5ad_to_seurat/
results/01.from_seurat_to_sce/
results/02.DropletUtils/
results/03.from_seurat_to_h5ad/
results/04.Scrublet/
results/05.scAutoQC/
results/06.Count_cells/
results/07.merge_cell_counts/
results/08.Generate_sample_report/
results/09.merge_dataset_reports/
```

Key final files:

| Output | Pattern |
|---|---|
| Final filtered h5ad | `results/05.scAutoQC/{dataset}/{dataset}_{sample}_scAutoQC_adata.h5ad` |
| Sample cell-count table | `results/06.Count_cells/{dataset}/{dataset}_{sample}_cell_counts.tsv` |
| Dataset cell-count summary | `results/07.merge_cell_counts/{dataset}_all_samples_cell_counts.tsv` |
| Sample QC report | `results/08.Generate_sample_report/{dataset}/{dataset}_{sample}_report.pdf` |
| Dataset merged report | `results/09.merge_dataset_reports/{dataset}_all_samples_report.pdf` |

---

## Citation

If you use scSIEVE-QC, please cite the underlying tools where appropriate:

- Snakemake
- DropletUtils
- Scrublet
- scAutoQC / sctk
- Seurat
- Scanpy

---

## License

This project is released under the MIT License. See `LICENSE` for details.
