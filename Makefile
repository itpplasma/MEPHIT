BUILD_DIR := build
BUILD_NINJA := $(BUILD_DIR)/build.ninja

# Prevent ambient shell variables from silently altering which libneo is fetched.
# unexport stops them reaching cmake via the child environment.
unexport LIBNEO_REF LIBNEO_PATH

# Forward only an explicitly passed command-line value, e.g. make LIBNEO_REF=main.
# $(origin ...) == "command line" only when the user typed it on the make invocation;
# "environment" (ambient shell) is deliberately excluded.
_LIBNEO_REF_FLAG :=
ifeq ($(origin LIBNEO_REF),command line)
  _LIBNEO_REF_FLAG := -DLIBNEO_REF=$(LIBNEO_REF)
endif

_LIBNEO_PATH_FLAG :=
ifeq ($(origin LIBNEO_PATH),command line)
  _LIBNEO_PATH_FLAG := -DLIBNEO_PATH=$(LIBNEO_PATH)
endif

.PHONY: all ninja test install clean
all: ninja

$(BUILD_NINJA):
	cmake --preset default -DCMAKE_COLOR_DIAGNOSTICS=ON $(_LIBNEO_REF_FLAG) $(_LIBNEO_PATH_FLAG)

ninja: $(BUILD_NINJA)
	cmake --build --preset default

test: ninja
	cd $(BUILD_DIR) && ctest

install: ninja
	cd $(BUILD_DIR) && ninja install

clean:
	rm -rf $(BUILD_DIR)
