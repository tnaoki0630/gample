import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ----- CSV 読み込み -----
df = pd.read_csv("sol.csv")
print(df)
# 必要な列を抽出
i = df["i"].values
j = df["j"].values

cols = ["rhs_decompose","sol_decompose","rhs_Ax","sol_Ax","rhs_Ay","sol_Ay"]

# 格子サイズを推定（i, j が 0 から開始している前提）
nx = i.max() + 1
ny = j.max() + 1

# 3 つのフィールドを (ny, nx) 配列に再構成
fields = {}
for c in cols:
    arr = np.full((ny, nx), np.nan)
    for _, row in df.iterrows():
        arr[int(row["j"]), int(row["i"])] = row[c]
    fields[c] = arr

# ----- ヒートマップを描画 -----
for c in cols:
    plt.figure(figsize=(6, 5))
    plt.title(c)
    plt.xlabel("i")
    plt.ylabel("j")
    plt.imshow(fields[c], origin="lower", aspect="auto")
    plt.colorbar(label=c)
    plt.show()
