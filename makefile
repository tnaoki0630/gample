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
TARGET = PIC
RESULT_DIR = bin

# 共通のソースファイル
SOURCES_COMMON = $(SRC_DIR)/main.mm \
                 $(SRC_DIR)/EMField.mm \
                 $(SRC_DIR)/Init.mm \
                 $(SRC_DIR)/Moment.mm \
                 $(SRC_DIR)/DebugPrint.mm\
                 $(SRC_DIR)/Constant.mm\

# EXEC の値によって Particle ソースを切り替え
ifeq ($(EXEC), mpi)
    PARTICLE_SRC = $(SRC_DIR)/Particle_mpi.mm
else
    PARTICLE_SRC = $(SRC_DIR)/Particle.mm
endif

# ソースリスト
SOURCES = $(SOURCES_COMMON) $(PARTICLE_SRC)

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
 
# デバッグ実行
debug: $(TARGET)
	chmod +x run_debug.sh
	./run_debug.sh $(TARGET) $(RESULT_DIR)

# clean
clean:
	rm -rf $(BUILD_DIR) $(TARGET)

.PHONY: all clean