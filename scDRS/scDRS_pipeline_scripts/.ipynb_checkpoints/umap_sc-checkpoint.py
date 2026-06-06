import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import warnings

warnings.filterwarnings("ignore")

# ==============================
# 基础设置
# ==============================

working_dir = '/biostack/home/suchenyi/graduate/scdrs'
os.chdir(working_dir)

sc.set_figure_params(dpi=150, facecolor='white')

os.makedirs("umap/S", exist_ok=True)
os.makedirs("umap/C", exist_ok=True)

print("Loading AnnData...")
adata = sc.read_h5ad("/biostack/home/suchenyi/scdrs/raw_hg.h5ad")

print("Discovering traits...")
score_files = sorted(glob.glob("trait_analysis/*.full_score.gz")) 
trait_list = [os.path.basename(f).replace(".full_score.gz", "") for f in score_files]

print(f"Found {len(trait_list)} traits.")

# ==============================
# Step 1: 读取所有 trait 到完整 adata
# ==============================

trait_mapping = {}
loaded_traits = []

for trait in trait_list:

    score_path = f"trait_analysis/{trait}.full_score.gz" 

    try:
        df = pd.read_csv(score_path, sep="\t", index_col=0)

        if "norm_score" not in df.columns:
            print(f"{trait}: no norm_score column")
            continue

        common_cells = adata.obs.index.intersection(df.index)

        if len(common_cells) == 0:
            print(f"{trait}: no overlapping cells")
            continue

        safe_trait = trait.replace(".", "_").replace("-", "_")
        adata.obs.loc[common_cells, safe_trait] = \
            df.loc[common_cells, "norm_score"].values

        trait_mapping[trait] = safe_trait
        loaded_traits.append(trait)

    except Exception as e:
        print(f"{trait} load error: {e}")
        continue

print(f"Successfully loaded {len(loaded_traits)} traits.")

# ==============================
# Step 2: 分区 
# ==============================

if "condition" not in adata.obs.columns:
    raise ValueError("adata.obs 中没有 'condition' 列")

adata_S = adata[adata.obs["condition"].str.contains("S", na=False)].copy()
adata_C = adata[adata.obs["condition"].str.contains("C", na=False)].copy()

print(f"S cells: {adata_S.n_obs}")
print(f"C cells: {adata_C.n_obs}")

# ==============================
# Step 3: 分区绘图
# ==============================

plot_count_S = 0
plot_count_C = 0

for trait in loaded_traits:

    col = trait_mapping[trait]

    for subset, label in [(adata_S, "S"), (adata_C, "C")]:

        if col not in subset.obs.columns:
            continue

        scores = subset.obs[col].dropna()

        if len(scores) == 0:
            continue

        vmin = np.percentile(scores, 1)
        vmax = np.percentile(scores, 99)

        if np.isnan(vmin) or np.isnan(vmax) or vmin == vmax:
            continue

        fig, ax = plt.subplots(figsize=(5, 4))

        try:
            sc.pl.embedding(
                subset,
                basis="umap.harmony",
                color=col,
                color_map="seismic",
                vcenter=0,
                vmin=vmin,
                vmax=vmax,
                size=2,
                ax=ax,
                show=False,
                frameon=False
            )

            ax.set_title(f"{trait}-GWAS ({label})", fontsize=10)

            safe_name = trait.replace("/", "_").replace("\\", "_")

            plt.savefig(
                f"umap/{label}/64human_{safe_name}_{label}.pdf",
                dpi=300,
                bbox_inches='tight'
            )
            plt.close()

            if label == "S":
                plot_count_S += 1
            else:
                plot_count_C += 1

        except Exception as e:
            print(f"{trait} ({label}) plot error: {e}")
            plt.close()
            continue

print("====================================")
print(f"S plots: {plot_count_S}")
print(f"C plots: {plot_count_C}")
print("Finished.")
