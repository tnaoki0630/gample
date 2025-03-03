#!/bin/bash

# 実行ファイルのパスを引数から取得
EXECUTABLE="$1"

# 実行ファイルのパスが空ならエラー
if [ -z "$EXECUTABLE" ]; then
    echo "Error: No executable provided."
    exit 1
fi

# 実行ファイルが存在しない場合のエラーチェック
if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable '$EXECUTABLE' not found."
    exit 1
fi

# 結果出力フォルダ内のデータを初期化
RESULT="$2"
if [[ "$RESULT_DIR" == *"/"* ]]; then
    echo "ディレクトリ名に '/' が含まれています。rm -r / が実行されるのを回避するため、エラー終了します。: $RESULT_DIR" >&2
    exit 1
fi
rm -r $RESULT
mkdir $RESULT

# Metal のデバッグ環境変数を設定
export METAL_DEVICE_WRAPPER_TYPE=1
export MTL_DEBUG_LAYER=1
export MTL_SHADER_VALIDATION=1

# 実行
echo "Running: $EXECUTABLE with Metal debug settings..."
./$EXECUTABLE