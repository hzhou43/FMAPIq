UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
    LIB_EXT = dylib
else ifeq ($(OS), Windows_NT)
    LIB_EXT = dll
else
    LIB_EXT = so
endif

SOURCES = ext/*.c
TARGET = src/sim2iq/libsim2iq.$(LIB_EXT)

all: $(TARGET)

$(TARGET): $(SOURCES)
	@mkdir -p src/sim2iq
	gcc -shared -fPIC -O3 -fopenmp -mfpmath=sse -Iext $(SOURCES) \
		-lfftw3f -lfftw3f_threads -lm -lgomp -o $(TARGET)
	@echo "✓ Built $(TARGET)"

build: $(TARGET)

clean:
	rm -f src/sim2iq/libsim2iq.*
	rm -f ext/*.o
	rm -rf build dist src/*.egg-info *.egg-info
	rm -f setup.cfg
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
	cd tests;./test.sh clean
	@echo "✓ Clean completed"

install: $(TARGET)
	pip3 install -e .

uninstall:
	pip3 uninstall -y sim2iq

reinstall: uninstall clean install

.PHONY: all build clean install uninstall reinstall
