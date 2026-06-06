import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import warnings
import glob

warnings.filterwarnings("ignore")

# 1. 设置工作路径和参数
sc.set_figure_params(dpi=150, facecolor='white')
working_dir = '/biostack/home/suchenyi/graduate/scdrs'
os.chdir(working_dir)

# 创建输出目录
os.makedirs("umap_hg38", exist_ok=True)

# 2. 读取原始数据
print("Loading AnnData...")
adata = sc.read_h5ad("/biostack/home/suchenyi/scdrs/raw_hg.h5ad")

# 3. 自动获取trait列表
print("Discovering trait results...")
score_files = glob.glob("score_hg38/*.full_score.gz")
trait_list = [os.path.basename(f).replace('.full_score.gz', '') for f in score_files]

if len(trait_list) == 0:
    print("Warning: No score files found in score_hg38/")
    print("Trying alternative pattern...")
    score_files = glob.glob("score_hg38/*.score.gz")
    trait_list = [os.path.basename(f).replace('.score.gz', '') for f in score_files]

print(f"Found {len(trait_list)} traits to process")

# 4. 读取 scDRS 结果并合并到 adata.obs
print("\nReading scDRS scores...")
conflict_traits = []
successful_traits = []
trait_mapping = {}  # 存储原始trait名到实际列名的映射

for trait in trait_list:  # 改为使用自动获取的trait_list
    score_path = f"score_hg38/{trait}.full_score.gz"
    
    # 添加文件存在性检查
    if not os.path.exists(score_path):
        print(f"Warning: File not found for {trait} at {score_path}")
        continue
        
    try:
        df_score = pd.read_csv(score_path, sep="\t", index_col=0)
        
        # 确保索引对齐
        common_cells = adata.obs.index.intersection(df_score.index)
        if len(common_cells) == 0:
            print(f"Error: No common cells found for {trait}. Check indices.")
            continue
            
        # 检查是否与基因名冲突
        if trait in adata.var_names:
            # 重命名，添加后缀以避免冲突
            new_trait_name = f"{trait}_score"
            conflict_traits.append((trait, new_trait_name))
            print(f"Warning: {trait} conflicts with gene name, renaming to {new_trait_name}")
            
            # 确保索引对齐后赋值
            adata.obs.loc[common_cells, new_trait_name] = df_score.loc[common_cells, "norm_score"].values
            trait_mapping[trait] = new_trait_name
        else:
            adata.obs.loc[common_cells, trait] = df_score.loc[common_cells, "norm_score"].values
            successful_traits.append(trait)
            trait_mapping[trait] = trait
            
    except Exception as e:
        print(f"Error loading {trait}: {e}")
        continue

print(f"\nSummary:")
print(f"Successfully loaded: {len(successful_traits)} traits")
print(f"Renamed due to conflicts: {len(conflict_traits)} traits")
for old_name, new_name in conflict_traits:
    print(f"  {old_name} -> {new_name}")

# 5. 绘制 dominant_labels (整体概览)
print("\nPlotting dominant labels...")
with plt.rc_context({"figure.figsize": (5, 5)}):
    sc.pl.embedding(
        adata,
        basis="umap.harmony",
        color="dominant_labels",
        show=False
    )
    plt.savefig("umap_hg38/64human_dominant_labels.pdf", bbox_inches='tight', dpi=300)
    plt.close()

# 6. 为每个 trait 单独绘制 UMAP
print("\nLooping through traits for plotting...")
plot_count = 0

for trait in trait_list:  # 改为使用自动获取的trait_list
    # 获取实际的列名（考虑可能的冲突重命名）
    col_name = trait_mapping.get(trait)
    if col_name is None or col_name not in adata.obs:
        print(f"Skipping {trait}: not loaded in adata.obs")
        continue
    
    # 检查是否有足够的有效值
    scores = adata.obs[col_name].dropna()
    if len(scores) == 0:
        print(f"Skipping {trait}: no valid scores")
        continue
    
    # 计算百分位数，处理极端值
    vmin = np.percentile(scores, 1)
    vmax = np.percentile(scores, 99)
    
    # 检查值范围是否有效
    if np.isnan(vmin) or np.isnan(vmax) or vmin == vmax:
        print(f"Skipping {trait}: invalid value range")
        continue
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(5, 4))
    
    try:
        # 使用vcenter使颜色映射对称
        sc.pl.embedding(
            adata,
            basis="umap.harmony",
            color=col_name,
            color_map="seismic",
            vcenter=0,  # 使0值在颜色映射中间
            vmin=vmin,
            vmax=vmax,
            size=2,
            ax=ax,
            show=False,
            frameon=False  # 去掉边框
        )
        
        ax.set_title(f"{trait}-GWAS", fontsize=10)
        
        # 保存图形
        output_name = trait.replace("/", "_").replace("\\", "_")  # 处理特殊字符
        plt.savefig(f"umap_hg38/64human_{output_name}_deepblue_bar.pdf", dpi=300, bbox_inches='tight')
        plt.close()
        
        plot_count += 1
        if plot_count % 10 == 0:
            print(f"Progress: {plot_count} traits plotted...")
            
    except Exception as e:
        print(f"Error plotting {trait}: {e}")
        plt.close()
        continue

print(f"\nAll plots completed! Total plots saved: {plot_count}")
print(f"Plots saved in 'umap_hg38/' folder.")
