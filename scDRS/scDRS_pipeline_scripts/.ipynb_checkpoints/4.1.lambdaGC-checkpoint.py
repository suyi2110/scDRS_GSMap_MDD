import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from scipy.stats import chi2

# =========================
# λGC函数
# =========================
def calc_lambda_gc(pvals):
    pvals = pvals.dropna()
    pvals = pvals[(pvals > 0)&(pvals < 1)]

    if len(pvals) == 0:
        return np.nan

    chisq = chi2.isf(pvals, 1)
    return np.median(chisq) / 0.4549


# =========================
# 参数
# =========================
input_dir = "/biostack/home/suchenyi/graduate/gwas/hg38"
output_prefix = "lambda_GC_all"

os.chdir(input_dir)

# =========================
# 自动获取所有tsv
# =========================
file_list = sorted(glob.glob("*.tsv"))

print(f"Found {len(file_list)} GWAS files")

results = []

print("========== λGC 计算 ==========")

for file in file_list:
    try:
        df = pd.read_csv(file, sep=r"\s+", engine="python")

        if 'P' not in df.columns:
            print(f"[跳过] {file}：没有 P 列")
            continue

        pvals = df['P']
        pvals = pd.to_numeric(pvals, errors='coerce')
        n_total = len(pvals)

        lambda_gc = calc_lambda_gc(pvals)
        lambda_gc_bg = calc_lambda_gc(pvals[pvals > 1e-5])

        trait = file.replace(".tsv", "")

        results.append({
            "Trait": trait,
            "N_SNP": n_total,
            "Lambda_GC": lambda_gc,
            "Lambda_GC_bg": lambda_gc_bg
        })

        print(f"{file:<25} λGC={lambda_gc:.4f}  bg={lambda_gc_bg:.4f}")

    except Exception as e:
        print(f"[错误] {file}: {e}")
        continue

print("=================================")

# =========================
# 保存结果
# =========================
df_res = pd.DataFrame(results)

df_res = df_res.sort_values("Lambda_GC", ascending=False)

df_res.to_csv(f"/biostack/home/suchenyi/graduate/scdrs/{output_prefix}.tsv", sep="\t", index=False)

print(f"结果已保存")

# =========================
# 画图（自动分页）
# =========================

def plot_bar(df, outname, top_n=30):
    plt.figure(figsize=(10,6))

    df_plot = df.head(top_n)

    x = np.arange(len(df_plot))

    plt.bar(x, df_plot["Lambda_GC"], label="λGC", alpha=0.8)
    plt.bar(x, df_plot["Lambda_GC_bg"], label="λGC_bg", alpha=0.7)

    plt.xticks(x, df_plot["Trait"], rotation=60, ha="right")
    plt.ylabel("λGC")
    plt.title(f"Top {top_n} GWAS λGC")
    plt.legend()

    plt.tight_layout()
    plt.savefig(outname, dpi=300)
    plt.close()

plot_bar(df_res, f"/biostack/home/suchenyi/graduate/scdrs/{output_prefix}.png")
print("绘图完成")
