#!/usr/bin/env python3
import json
import copy
from pathlib import Path

# テンプレート JSON と出力ディレクトリ
input_file = Path("condition.json")

# テンプレート読み込み
with input_file.open(encoding="utf-8") as f:
    template = json.load(f)

# 0～8 のループ
for i in range(9):
    for j in range(i+1):
        tgs = 2**j
        ics = 2**i
        proj = f"tgs{tgs}ics{ics}"

        # テンプレートを深いコピーして書き換え
        cfg = copy.deepcopy(template)
        cfg["ParamForTimeIntegration"]["ProjectName"]    = proj
        cfg["ParamForComputing"]["ThreadGroupSize"]       = tgs
        cfg["ParamForComputing"]["IntegrationChunkSize"]  = ics

        # ファイルに書き出し
        out_path = Path(f"{proj}.json")
        with out_path.open("w", encoding="utf-8") as wf:
            json.dump(cfg, wf, indent=2, ensure_ascii=False)

        print(f"Generated: {out_path}")

