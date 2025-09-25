#!/bin/bash
# Build executable for FMAPIq Pipeline

# Pin Python version for consistent builds
uv python pin 3.8

# Add GUI dependency
uv add gradio==3.50.2

# Build executable with PyInstaller
uv run pyinstaller \
    --onefile \
    --clean \
    --name fmapiq \
    main.py \
    --collect-data gradio_client \
    --collect-data gradio \
    --collect-data fmapb2 \
    --collect-data maeppi \
    --collect-data sim2iq \
    --add-binary .venv/lib/python3.8/site-packages/fmapb2/*.so:fmapb2 \
    --add-binary .venv/lib/python3.8/site-packages/maeppi/*.so:maeppi \
    --add-binary .venv/lib/python3.8/site-packages/sim2iq/*.so:sim2iq