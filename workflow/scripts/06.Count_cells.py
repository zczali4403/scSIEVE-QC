import sys
import scanpy as sc
import pandas as pd
from pathlib import Path

# 输入参数：样本名、dropletutils文件、scrublet文件、scAutoQC文件、输出结果文件
sample = sys.argv[1]
droplet_file = sys.argv[2]
scrublet_file = sys.argv[3]
scautoqc_file = sys.argv[4]
outfile = sys.argv[5]

# 统计 cell 数
counts = {}
counts["sample"] = sample

# droplet
adata = sc.read_h5ad(droplet_file)
counts["after_dropletutils"] = adata.n_obs

# scrublet
adata = sc.read_h5ad(scrublet_file)
counts["after_scrublet"] = adata.n_obs

# scAutoQC
adata = sc.read_h5ad(scautoqc_file)
counts["after_scautoqc"] = adata.n_obs

# 保存
df = pd.DataFrame([counts])
df.to_csv(outfile, sep="\t", index=False)
