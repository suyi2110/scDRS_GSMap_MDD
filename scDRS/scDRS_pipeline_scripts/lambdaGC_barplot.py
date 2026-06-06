import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# =========================
# 读取数据
# =========================
file = "/biostack/home/suchenyi/graduate/scdrs/lambda_GC_all.tsv"
df = pd.read_csv(file, sep="\t")

# 去掉NaN
df = df.dropna(subset=["Lambda_GC"])

# 按λGC排序
df = df.sort_values("Lambda_GC", ascending=False)
# =========================
# 图1：全部traits柱状图（竖版）
# =========================
plt.figure(figsize=(18,6))

x = np.arange(len(df))

plt.bar(x, df["Lambda_GC"], label="λGC", alpha=0.8)
plt.bar(x, df["Lambda_GC_bg"], label="λGC_bg", alpha=0.7)

plt.xticks(x, df["Trait"], rotation=90, fontsize=6)
plt.ylabel("λGC")
plt.title("All GWAS Traits λGC")

plt.legend()
plt.tight_layout()

plt.savefig("/biostack/home/suchenyi/graduate/scdrs/lambdaGC_all_bar.png", dpi=300)
plt.close()

print("全量柱状图已保存")
# =========================
#横向
# =========================
plt.figure(figsize=(8, max(6, len(df)*0.25)))

y = np.arange(len(df))

plt.barh(y, df["Lambda_GC"], label="λGC", alpha=0.8)
plt.barh(y, df["Lambda_GC_bg"], label="λGC_bg", alpha=0.7)

plt.yticks(y, df["Trait"], fontsize=6)
plt.xlabel("λGC")
plt.title("All GWAS Traits λGC (sorted)")

plt.legend()
plt.tight_layout()

plt.savefig("/biostack/home/suchenyi/graduate/scdrs/lambdaGC_all_barh.png", dpi=300)
plt.close()

print("横向柱状图已保存")