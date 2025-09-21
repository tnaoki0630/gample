import numpy as np
import re
import glob
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

# バイナリファイルのレコード形式（time:int, x, y, vx, vy, vz:double）
record_dtype = np.dtype([
    ("time", np.int32),
    ("x", np.float32),
    ("y", np.float32),
    ("vx", np.float32),
    ("vy", np.float32),
    ("vz", np.float32)
])

args = sys.argv
## projectname
pname = args[1]

# ファイル名のパターンに合わせて全粒子のファイルを取得
files_unsorted = glob.glob(f"bin/{pname}_PhaseSpace_electron_*.bin")
files = sorted(files_unsorted, key=lambda s: int(re.search(r'\d+', s).group()))
print(files)
# files = sorted(glob.glob("bin/large_PhaseSpace_ion_Xe1_*.bin"))
if not files:
    raise FileNotFoundError("*.bin が見つかりません。")

# 各粒子のデータを読み込み（各ファイルは1粒子分の時系列データ）
particle_data = [np.fromfile(f, dtype=record_dtype) for f in files]
print(particle_data)

num_particles = len(particle_data)
num_frames = particle_data[0].shape[0]
for d in particle_data:
    if d.shape[0] != num_frames:
        raise ValueError("粒子間でフレーム数が一致していません。")

# 各粒子の (x, y) をベクトル化して frame_positions を作成（shape: (num_frames, num_particles, 2)）
frame_positions = np.stack([np.column_stack((pd["x"], pd["y"])) for pd in particle_data], axis=1)
# 時刻情報は最初の粒子のデータから取得
times = particle_data[0]["time"]

# プロットの設定
fig, ax = plt.subplots()
scat = ax.scatter(frame_positions[0, :, 0], frame_positions[0, :, 1], s=5)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title(f"Time step: {times[0]}")
ax.set_xlim(0, 0.1)
ax.set_ylim(0, 0.1)

# ax.axvline(2.4, color='gray', linestyle=':', linewidth=0.5) ## cathodeline
ax.set_aspect('equal', adjustable='box')

# 粒子番号をプロット
x = frame_positions[0, :, 0]
y = frame_positions[0, :, 1]
for i, (xx, yy) in enumerate(zip(x, y)):
    print(i,xx,yy)
    ax.text(xx, yy+0.002, f"{i}", ha='center', va='center', fontsize=6, clip_on=True)
plt.show()

def update(frame):
    # フレームごとに粒子の位置を更新
    scat.set_offsets(frame_positions[frame])
    ax.set_title(f"Time step: {times[frame]}")
    return scat,

# アニメーション作成
ani = animation.FuncAnimation(fig, update, frames=num_frames, interval=100, blit=True)

# 動画として保存
ani.save(f"fig/{pname}_ptraj.mp4", writer="ffmpeg", fps=60)