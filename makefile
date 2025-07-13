# require `brew install libomp`
OMP_INC = $(shell brew --prefix libomp)/include
OMP_LIB = $(shell brew --prefix libomp)/lib
BOOST_INC = $(shell brew --prefix boost)/include

# コンパイラとフラグの設定
CXX = clang++
CXXFLAGS = -std=c++17 -O3\
           -x objective-c++ \
           -isysroot $(shell xcrun --show-sdk-path)\
           -I/Users/$(shell whoami)/Library/amgcl\
           -I/opt/homebrew/Cellar/boost/1.87.0/include\
           -Xpreprocessor -fopenmp \
           -I$(OMP_INC) -I$(BOOST_INC)

# リンカフラグ
LDFLAGS = -framework Metal \
          -framework Foundation \
          -framework MetalKit \
          -fobjc-arc\
          -L$(shell xcrun --show-sdk-path)/usr/lib\
          -F$(shell xcrun --show-sdk-path)/System/Library/Frameworks\
          -lomp -L$(OMP_LIB)

# ディレクトリ設定
SRC_DIR = src
BUILD_DIR = build
TARGET = PIC
RESULT_DIR = bin

# 共通のソースファイル
SOURCES_COMMON = $(SRC_DIR)/main.mm \
                 $(SRC_DIR)/Init.mm \
                 $(SRC_DIR)/EMField.mm \
                 $(SRC_DIR)/Moment.mm \
                 $(SRC_DIR)/FileIO.mm \
                 $(SRC_DIR)/Collision.mm \

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

# clean
clean:
	rm -rf $(BUILD_DIR) $(TARGET)

.PHONY: all clean
