BUILD_DIR := build
MAKEFILE := $(BUILD_DIR)/Makefile
NUM_CORES := $(shell nproc)

.PHONY: all build test install clean
all: build

$(MAKEFILE):
	cmake --preset default

build: $(MAKEFILE)
	cmake --build --preset default --parallel $(NUM_CORES)

test: build
	cd $(BUILD_DIR) && ctest

install: build
	cd $(BUILD_DIR) && make install

clean:
	rm -rf $(BUILD_DIR)
