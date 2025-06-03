import matplotlib.pyplot as plt

# 読み込むファイル名
filename = 'stdout.txt'

# ファイルを開いて各行の数値をリストに格納
values = []
with open(filename, 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        try:
            x = float(line)
            values.append(x)
        except ValueError:
            # 数値に変換できない行はスキップ
            continue

# ヒストグラムを作成
plt.figure(figsize=(6, 4))
plt.hist(values, bins=100, range=(0, 2.5), edgecolor='black')
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Histogram (0 ～ 2.5)')

# グラフを表示
plt.tight_layout()
plt.show()
