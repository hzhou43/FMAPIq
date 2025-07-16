UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
    LIB_EXT = dylib
else ifeq ($(OS), Windows_NT)
    LIB_EXT = dll
else
    LIB_EXT = so
endif

SOURCES = ext/fftw.c ext/grd.c ext/map.c ext/score.c
TARGET = src/fmapb2/libfmapb2.$(LIB_EXT)

all: $(TARGET)

$(TARGET): $(SOURCES)
	@mkdir -p src/fmapb2
	gcc -shared -fPIC -O3 -fopenmp -mfpmath=sse -Iext $(SOURCES) \
		-lfftw3 -lfftw3_threads -lm -lgomp -o $(TARGET)
	@echo "✓ Built $(TARGET)"

build: $(TARGET)

clean:
	rm -f src/fmapb2/libfmapb2.*
	rm -f ext/*.o
	rm -rf build dist src/*.egg-info *.egg-info
	rm -f setup.cfg
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
	cd tests;./test.sh clean
	@echo "✓ Clean completed"

install: $(TARGET)
	pip3 install .

dev-install: $(TARGET)
	pip install -e .

uninstall:
	pip3 uninstall -y fmapb2

reinstall: uninstall clean install

.PHONY: all build clean install uninstall reinstall
