BUILD_DIR := build
BUILD_NINJA := $(BUILD_DIR)/build.ninja

.PHONY: all ninja test install clean
all: ninja

$(BUILD_NINJA):
	cmake --preset default -DCMAKE_COLOR_DIAGNOSTICS=ON

ninja: $(BUILD_NINJA)
	cmake --build --preset default

test: ninja
	cd $(BUILD_DIR) && ctest

install: ninja
	cd $(BUILD_DIR) && ninja install

clean:
	rm -rf $(BUILD_DIR)
