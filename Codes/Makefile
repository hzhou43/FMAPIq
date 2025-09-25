.PHONY: all fmapb2 maeppi sim2iq clean

# Default target - build all packages
all: fmapb2 maeppi sim2iq

# Build individual packages
fmapb2:
	$(MAKE) -C fmapb2-pkg

maeppi:
	$(MAKE) -C maeppi-pkg

sim2iq:
	$(MAKE) -C sim2iq-pkg

# Clean all packages
clean:
	$(MAKE) -C fmapb2-pkg clean
	$(MAKE) -C maeppi-pkg clean
	$(MAKE) -C sim2iq-pkg clean
