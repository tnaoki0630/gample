# Makefile for integrated Metal reduction sample (Objective-C++)

APP = main
SRCS = main.mm
OBJS = $(SRCS:.mm=.o)
CC = clang++
CXXFLAGS = -std=c++17 \
           -isysroot $(shell xcrun --show-sdk-path) \
           -I/Users/$(shell whoami)/Library/amgcl\
           -I/opt/homebrew/Cellar/boost/1.87.0/include


# リンカフラグ
LDFLAGS = -framework Metal \
          -framework Foundation \
          -framework MetalKit \
          -fobjc-arc\
          -L$(shell xcrun --show-sdk-path)/usr/lib\
          -F$(shell xcrun --show-sdk-path)/System/Library/Frameworks

all: $(APP)

$(APP): $(OBJS)
	$(CC) $(CXXFLAGS) -o $(APP) $(OBJS) $(LDFLAGS)

%.o: %.mm
	$(CC) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(APP)

.PHONY: all clean
