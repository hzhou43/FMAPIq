#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FMAPIq Pipeline - Unified Launcher

Usage:
    python main.py                    # Launch GUI (if available)
    python main.py subA.pqr [opts]    # CLI mode
    python main.py --gui              # Force GUI mode
    python main.py --cli              # Force CLI mode
"""

import sys

def main():
    """Main entry point that decides between GUI and CLI modes"""
    
    # Check for explicit mode flags
    force_gui = '--gui' in sys.argv
    force_cli = '--cli' in sys.argv
    
    # Remove mode flags from arguments
    if force_gui:
        sys.argv.remove('--gui')
    if force_cli:
        sys.argv.remove('--cli')
    
    # Determine execution mode
    if force_gui or (len(sys.argv) == 1 and not force_cli):
        launch_gui()
    else:
        launch_cli()

def launch_gui():
    """Launch the GUI interface"""
    try:
        import gradio
        from gui import main as gui_main
        print("Launching GUI interface...")
        gui_main()
    except ImportError:
        print("GUI not available. Install gradio: pip install gradio")
        print("Falling back to CLI mode...")
        print(__doc__.strip())
        show_cli_examples()

def launch_cli():
    """Launch the command line interface"""
    from pipeline import main as cli_main
    cli_main()

def show_cli_examples():
    """Display CLI usage examples"""
    print("\nCLI Examples:")
    print("  python main.py subA.pqr                    # Run full pipeline")
    print("  python main.py --load /path/to/workdir     # Load existing work")
    print("  python main.py subA.pqr --pack             # Run and pack results")
    print("  python main.py --help                      # Show all options")

if __name__ == "__main__":
    main()
