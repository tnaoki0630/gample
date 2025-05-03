#!/bin/bash

# 既存のバイナリファイルを削除
rm ./bin/*.bin
mkdir ./bin

# Metal のデバッグ環境変数を設定
export METAL_DEVICE_WRAPPER_TYPE=1
export MTL_DEBUG_LAYER=1
export MTL_SHADER_VALIDATION=1

# 実行
echo "Running: PIC with Metal debug settings..."
./PIC