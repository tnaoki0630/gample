# コンパイラとフラグの設定
CXX = clang++
CXXFLAGS = -std=c++17 \
           -isysroot $(shell xcrun --show-sdk-path) 

# リンカフラグ
LDFLAGS = -framework Metal \
          -framework Foundation \
          -framework MetalKit \
          -fobjc-arc\
          -L$(shell xcrun --show-sdk-path)/usr/lib\
          -F$(shell xcrun --show-sdk-path)/System/Library/Frameworks

# ディレクトリ設定
SRC_DIR = src
BUILD_DIR = build
TARGET = PIC.out

# ソースファイル
SOURCES = $(SRC_DIR)/PIC.mm \
          $(SRC_DIR)/main.mm

# オブジェクトファイル
OBJECTS = $(SOURCES:$(SRC_DIR)/%.mm=$(BUILD_DIR)/%.o)

# デフォルトターゲット
all: $(BUILD_DIR) $(TARGET)

# ビルドディレクトリの作成
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# 実行ファイルの生成
$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@

# オブジェクトファイルの生成
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.mm
	$(CXX) $(CXXFLAGS) -c $< -o $@

# クリーン
clean:
	rm -rf $(BUILD_DIR) $(TARGET)

.PHONY: all clean