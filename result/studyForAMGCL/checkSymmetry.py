import scipy.sparse as sp
from scipy.io import mmread
import numpy as np

# 行列を読み込む
A = mmread("poisson.mtx").tocoo()

# 対称性チェック
diff = (A - A.T).tocoo()
max_diff = np.abs(diff.data).max() if diff.data.size > 0 else 0.0
print("非対称成分の最大絶対値 =", max_diff)

if max_diff < 1e-12:
    print("→ 十分小さいので、実質的に対称行列とみなせます。")
else:
    print("→ 要注意：差分が {:.3e} なので非対称成分が残っています。".format(max_diff))
