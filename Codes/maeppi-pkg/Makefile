UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
    LIB_EXT = dylib
else ifeq ($(OS), Windows_NT)
    LIB_EXT = dll
else
    LIB_EXT = so
endif

FMAP_SRC = ext/fmap_energy.c
SIM_SRC = ext/simulation.c ext/particle_simple.c ext/fmap_energy.c
FMAP_TARGET = src/maeppi/libmaeppiern.$(LIB_EXT)
SIM_TARGET = src/maeppi/libmaeppisim.$(LIB_EXT)

all: $(FMAP_TARGET) $(SIM_TARGET)
#-mfpmath=sse
$(FMAP_TARGET): $(FMAP_SRC)
	gcc -shared -fPIC -O3 -fopenmp -Iext $(FMAP_SRC) -lm -lgomp -o $(FMAP_TARGET)
	@echo "✓ Built $(FMAP_TARGET)"

$(SIM_TARGET): $(SIM_SRC)
	gcc -shared -fPIC -O3 -fopenmp -Iext $(SIM_SRC) -lm -lgomp -o $(SIM_TARGET)
	@echo "✓ Built $(SIM_TARGET)"

build: $(FMAP_TARGET) $(SIM_TARGET)

clean:
	rm -f src/maeppi/libmaeppiern.*
	rm -f src/maeppi/libmaeppisim.*
	rm -f ext/*.o
	rm -rf build dist src/*.egg-info *.egg-info
	rm -f setup.cfg
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
	cd tests;./test.sh clean;./run.fmapb2.sh clean
	@echo "✓ Clean completed"

install: $(FMAP_TARGET) $(SIM_TARGET)
	pip3 install .

dev-install: $(TARGET)
	pip3 install -e .

uninstall:
	pip3 uninstall -y maeppi

reinstall: uninstall clean install

.PHONY: all build clean install uninstall reinstall
