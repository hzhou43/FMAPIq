#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FMAPIq Pipeline - Gradio GUI Wrapper with Session Management
"""

import os
import sys
import glob
import uuid
import tempfile

try:
    import gradio as gr
except ImportError:
    print("Gradio not available. Install with: pip install gradio")
    sys.exit(1)

try:
    from pipeline import Pipeline, DEFAULTS
except ImportError:
    print("Error: pipeline.py not found. Make sure it's in the same directory.")
    sys.exit(1)

class SessionManager:
    """Manages multiple pipeline instances with unique work directories"""
    
    def __init__(self, base_dir="sessions"):
        self.base_dir = base_dir
        self.sessions = {}
        self._create_base_dir()
    
    def _create_base_dir(self):
        """Create base sessions directory if it doesn't exist"""
        if not os.path.exists(self.base_dir):
            os.makedirs(self.base_dir)
    
    def create_session(self, session_id=None):
        """Create new session with unique ID and work directory"""
        if not session_id:
            session_id = str(uuid.uuid4())[:8]
        
        if session_id in self.sessions:
            return session_id, self.sessions[session_id]
        
        work_dir = os.path.join(self.base_dir, session_id)
        pipeline = Pipeline(work_dir)
        self.sessions[session_id] = pipeline
        
        return session_id, pipeline
    
    def get_session(self, session_id):
        """Get existing session by ID"""
        return self.sessions.get(session_id)
    
    def list_sessions(self):
        """List all active session IDs"""
        return list(self.sessions.keys())
    
    def remove_session(self, session_id):
        """Remove session from memory (files are preserved)"""
        if session_id in self.sessions:
            del self.sessions[session_id]

def create_gui():
    """Create and configure the Gradio interface"""
    
    session_manager = SessionManager()
    current_session = {"id": None, "pipeline": None}

    def ensure_session():
        """Auto-create session if none exists"""
        if not current_session["pipeline"]:
            session_id, pipeline = session_manager.create_session()
            current_session["id"] = session_id
            current_session["pipeline"] = pipeline
            return f"Auto-created session: {session_id}\n"
        return ""

    def get_current_pipeline():
        """Get current pipeline, creating session if needed"""
        ensure_session()
        return current_session["pipeline"]

    def get_session_display():
        """Get current session ID for display"""
        if current_session["id"]:
            return f'Current Session ID: {current_session["id"]}'
        return "Current Session ID: None"

    def calculate_concentrations(n_proteins, box_size):
        """Calculate molecular weight and concentrations"""
        pipeline = get_current_pipeline()
        if not pipeline.work_dir or not hasattr(pipeline, 'mass') or pipeline.mass is None:
            return "N/A", "N/A", "N/A"
    
        try:
            # Convert box size from Angstroms to cm
            box_size_cm = box_size * 1e-8
            volume_cm3 = box_size_cm ** 3
    
            mw_da = pipeline.mass
    
            # Calculate molar concentration
            avogadro = 6.022e23
            moles = n_proteins / avogadro
            volume_liters = volume_cm3 * 1e-3
    
            molarity = moles / volume_liters
            concentration_mm = molarity * 1000
    
            # Calculate mass concentration
            mass_grams = moles * mw_da
            concentration_mg_ml = (mass_grams * 1000) / volume_cm3
    
            return f"{mw_da:.1f}", f"{concentration_mg_ml:.2f}", f"{concentration_mm:.3f}"
    
        except Exception as e:
            print(f"Error in concentration calculation: {e}")
            return "Error", "Error", "Error"

    def update_concentrations(n_proteins, box_size):
        """Update concentration displays when parameters change"""
        return calculate_concentrations(n_proteins, box_size)

    def load_session_or_archive(session_id, archive_file):
        """Load session or unpack archive file"""
        auto_msg = ensure_session()
        
        # Handle archive file upload
        if archive_file is not None:
            filename = archive_file.name.lower()
            if not (filename.endswith('.tgz') or filename.endswith('.tar.gz') or filename.endswith('.gz')):
                return f"{auto_msg}Error: Invalid file type. Only .tgz and .tar.gz files are supported.\nReceived: {archive_file.name}", get_session_display()
            
            # Create new pipeline for archive
            new_session_id, pipeline = session_manager.create_session()
            work_dir, message = pipeline._unpack_tgz(archive_file.name)
            if not work_dir:
                return f"{auto_msg}{message}", get_session_display()
            
            pipeline.work_dir = work_dir
            current_session["id"] = new_session_id
            current_session["pipeline"] = pipeline
            status = f"{auto_msg}Unpacked archive: {os.path.basename(archive_file.name)}\nSession: {new_session_id}\nWork directory: {work_dir}\n\n"
            
        # Handle session ID
        elif session_id and session_id.strip():
            pipeline = session_manager.get_session(session_id)
            if not pipeline:
                session_id, pipeline = session_manager.create_session(session_id)
                status = f"Created new session: {session_id}"
            else:
                status = f"Loaded existing session: {session_id}"
            
            current_session["id"] = session_id
            current_session["pipeline"] = pipeline
            
            if pipeline.work_dir and os.path.exists(pipeline.work_dir):
                status += f"\nWork directory: {pipeline.work_dir}"
            else:
                status += f"\nWork directory: {pipeline.work_dir} (will be created)"
                
        else:
            return f"{auto_msg}Please either select a session ID or upload an archive file", get_session_display()

        # Load pipeline data if available
        if pipeline.work_dir and os.path.exists(pipeline.work_dir):
            if pipeline.load_data():
                stage1_complete = os.path.exists(os.path.join(pipeline.work_dir, "parms.txt"))
                stage2_complete = len(glob.glob(os.path.join(pipeline.work_dir, "config_*.dat"))) > 0
                stage3_complete = os.path.exists(os.path.join(pipeline.work_dir, "Iqs.txt"))

                status += f"\nStage 1: {'OK' if stage1_complete else 'X'}"
                status += f"\nStage 2: {'OK' if stage2_complete else 'X'}"
                status += f"\nStage 3: {'OK' if stage3_complete else 'X'}"

                if stage3_complete and pipeline.iq_data:
                    status += f"\nLoaded {pipeline.iq_data['n_curves']} I(q) datasets"
            else:
                status += "\nFailed to load pipeline data"

        return status, get_session_display()

    def run_stage1(pqr_file, ion, temp, escl, vscl, threads_val, update_only):
        """Run or update Stage 1: Energy Maps"""
        auto_msg = ensure_session()
        pipeline = current_session["pipeline"]
        
        if update_only:
            if pipeline.load_data(stage=1):
                return f"{auto_msg}Stage 1 Updated!\n{pipeline.params}", get_session_display()
            else:
                return f"{auto_msg}Failed to load Stage 1 data", get_session_display()
        elif pqr_file:
            result = pipeline.stage1(pqr_file.name, ion, temp, escl, vscl, int(threads_val))
            return f"{auto_msg}{result}", get_session_display()
        elif pipeline.pqr_file:
            result = pipeline.stage1(None, ion, temp, escl, vscl, int(threads_val))
            return f"{auto_msg}{result}", get_session_display()
        return f"{auto_msg}No PQR file provided", get_session_display()

    def run_stage2(n_reps, n_cycles, n_proteins, box_size, threads, sim_threads, update_only):
        """Run or update Stage 2: Simulation"""
        auto_msg = ensure_session()
        pipeline = current_session["pipeline"]
        
        if update_only:
            if pipeline.load_data(stage=2):
                return f"{auto_msg}Stage 2 Updated! Found {len(pipeline.configs)} configurations", get_session_display()
            else:
                return f"{auto_msg}Failed to load Stage 2 data", get_session_display()
        else:
            result = pipeline.stage2(int(n_reps), int(n_cycles), int(n_proteins),
                                 box_size, int(threads), int(sim_threads))
            return f"{auto_msg}{result}", get_session_display()

    def run_stage3(box_size, spacing, threads, update_only):
        """Run or update Stage 3: Calculate I(q)"""
        auto_msg = ensure_session()
        pipeline = current_session["pipeline"]
        
        if update_only:
            if pipeline.load_data(stage=3):
                if pipeline.iq_data:
                    return f"{auto_msg}Stage 3 Updated! Loaded {pipeline.iq_data['n_curves']} I(q) replicates", get_session_display()
                else:
                    return f"{auto_msg}Stage 3 Updated! No I(q) data found.", get_session_display()
            else:
                return f"{auto_msg}Failed to load Stage 3 data", get_session_display()
        else:
            result = pipeline.stage3(box_size, spacing, threads)
            return f"{auto_msg}{result}", get_session_display()

    def create_plot(show_individual, show_pq, q_min, q_max, y_min, y_max, log_scale_x, log_scale_y):
        """Create I(q) plot"""
        auto_msg = ensure_session()
        pipeline = current_session["pipeline"]
        plot_path, msg = pipeline.create_plot(show_individual, show_pq,
                                            q_min, q_max, y_min, y_max, log_scale_x, log_scale_y)
        return plot_path

    def pack_results(skip_ve, skip_raw):
        """Pack work directory for download"""
        auto_msg = ensure_session()
        pipeline = current_session["pipeline"]
        file_path, message = pipeline.pack_results(skip_ve, skip_raw)
        return file_path, f"{auto_msg}{message}"

    def get_stats_file():
        """Get path to statistics file if it exists"""
        pipeline = get_current_pipeline()
        if pipeline.work_dir:
            stats_path = os.path.join(pipeline.work_dir, "Iqs.stat.txt")
            if os.path.exists(stats_path):
                return stats_path
        return None

    def get_iqs_file():
        """Get path to I(q) data file if it exists"""
        pipeline = get_current_pipeline()
        if pipeline.work_dir:
            iqs_path = os.path.join(pipeline.work_dir, "Iqs.txt")
            if os.path.exists(iqs_path):
                return iqs_path
        return None

    def update_jobs_display(total_threads, sim_threads_val):
        """Calculate and display number of parallel jobs"""
        jobs = max(1, int(total_threads) // int(sim_threads_val))
        return f"{jobs} (= {int(total_threads)} / {int(sim_threads_val)})"

    def toggle_load_section(show):
        """Show/hide session management section"""
        return gr.Group(visible=show)

    # Create the Gradio interface
    with gr.Blocks(title="FMAPIq Pipeline") as app:
        gr.Markdown("## FMAPIq Pipeline")

        # Session management section
        with gr.Row():
            show_load_option = gr.Checkbox(label="Show session management", value=False)
            session_display = gr.Textbox(
                label="Current Session ID", 
                value="Current Session ID: None",
                interactive=False,
                container=False
            )

        with gr.Group(visible=False) as load_section:
            with gr.Row():
                with gr.Column():
                    session_dropdown = gr.Dropdown(
                        choices=session_manager.list_sessions(),
                        label="Load Existing Session ID",
                        value=None,
                        allow_custom_value=True
                    )
                    archive_upload = gr.File(
                        label="Or Load from Archive",
                        file_types=[".tgz", ".tar.gz", ".gz"],
                        file_count="single"
                    )
                with gr.Column():
                    load_btn = gr.Button("Load")
                    load_status = gr.Textbox(label="Status", lines=4)

        # Stage 1: Energy Maps
        with gr.Tab("Stage 1: Energy Maps"):
            with gr.Row():
                with gr.Column():
                    pqr_file = gr.File(label="PQR File", file_types=[".pqr"])
                    with gr.Row():
                        ion_strength = gr.Number(value=DEFAULTS['ion'], label="Ion Strength (M)")
                        temperature = gr.Number(value=DEFAULTS['temp'], label="Temperature (C)")
                    with gr.Row():
                        escl = gr.Number(value=DEFAULTS['escl'], label="escl")
                        vscl = gr.Number(value=DEFAULTS['vscl'], label="vscl (Negative for default based on Mw)")
                    threads = gr.Slider(1, os.cpu_count(), value=os.cpu_count()//2, step=1, label="Threads")
                    update_only_1 = gr.Checkbox(label="Update Only", value=False)
                    stage1_btn = gr.Button("Run Stage 1")
                with gr.Column():
                    stage1_output = gr.Textbox(label="Results", lines=8)

        # Stage 2: Simulation
        with gr.Tab("Stage 2: Simulation"):
            with gr.Row():
                with gr.Column():
                    gr.Markdown("### Concentration options")
                    with gr.Row():
                        n_proteins = gr.Number(value=DEFAULTS['n_proteins'], label="Number of proteins", precision=0)
                        box_size = gr.Number(value=DEFAULTS['box_size'], label="Box Size (Å)")
                    
                    # Concentration display
                    with gr.Row():
                        mw_display = gr.Textbox(
                            value="N/A", 
                            label="Molecular Weight (Da)", 
                            interactive=False
                        )
                        conc_mg_ml_display = gr.Textbox(
                            value="N/A", 
                            label="Concentration (mg/mL)", 
                            interactive=False
                        )
                        conc_mm_display = gr.Textbox(
                            value="N/A", 
                            label="Concentration (mM)", 
                            interactive=False
                        )
                    
                    gr.Markdown("### Simulation options")
                    with gr.Row():
                        n_reps = gr.Number(value=DEFAULTS['n_reps'], label="Replicates", precision=0)
                        n_cycles = gr.Number(value=DEFAULTS['n_cycles'], label="Cycles", precision=0)
                    
                    with gr.Row():
                        sim_threads = gr.Slider(1, os.cpu_count()*2, value=1, step=1, label="Threads per Job")
                        jobs_display = gr.Textbox(
                            value=f"{max(1, os.cpu_count()//2 // 1)}",
                            label="Parallel Jobs",
                            interactive=False
                        )
                    update_only_2 = gr.Checkbox(label="Update Only", value=False)
                    stage2_btn = gr.Button("Run Stage 2")
                with gr.Column():
                    stage2_output = gr.Textbox(label="Results", lines=8)

        # Stage 3: Calculate I(q)
        with gr.Tab("Stage 3: Calculate I(q)"):
            with gr.Row():
                with gr.Column():
                    spacing = gr.Number(value=DEFAULTS['spacing'], label="Spacing")
                    update_only_3 = gr.Checkbox(label="Update Only", value=False)
                    stage3_btn = gr.Button("Calculate I(q)")
                with gr.Column():
                    stage3_output = gr.Textbox(label="Results", lines=8)

        # Plot I(q)
        with gr.Tab("Plot I(q)"):
            with gr.Row():
                with gr.Column():
                    with gr.Row():
                        show_individual = gr.Checkbox(label="Show individual curves", value=DEFAULTS['show_individual'])
                        show_pq = gr.Checkbox(label="Show P(q) reference", value=DEFAULTS['show_pq'])
                    with gr.Row():
                        log_scale_x = gr.Checkbox(label="Log scale (x-axis)", value=DEFAULTS['log_scale_x'])
                        log_scale_y = gr.Checkbox(label="Log scale (y-axis)", value=DEFAULTS['log_scale_y'])
                    with gr.Row():
                        q_min = gr.Number(value=DEFAULTS['q_min'], label="Q min")
                        q_max = gr.Number(value=DEFAULTS['q_max'], label="Q max")
                    with gr.Row():
                        y_min = gr.Number(value=DEFAULTS['y_min'], label="Y min")
                        y_max = gr.Number(value=DEFAULTS['y_max'], label="Y max")
                    plot_btn = gr.Button("Generate Plot")
                with gr.Column():
                    plot_image = gr.Image(label="I(q) Plot")

        # Download
        with gr.Tab("Download"):
            with gr.Row():
                with gr.Column():
                    gr.Markdown("### Direct File Downloads")
                    stats_file = gr.File(label="Download Iqs.stat.txt", interactive=False)
                    iqs_file = gr.File(label="Download Iqs.txt", interactive=False)
                with gr.Column():
                    gr.Markdown("### Pack Work Directory")
                    skip_ve = gr.Checkbox(label="Skip energy maps", value=True)
                    skip_raw = gr.Checkbox(label="Skip raw data", value=True)
                    pack_btn = gr.Button("Pack Work Directory")
                    pack_status = gr.Textbox(label="Pack Status", lines=2)
                    download_file = gr.File(label="Download Archive", interactive=False)

        # Event handlers
        show_load_option.change(toggle_load_section, [show_load_option], [load_section])
        
        load_btn.click(
            load_session_or_archive, 
            [session_dropdown, archive_upload], 
            [load_status, session_display]
        )
        
        stage1_btn.click(run_stage1,
                      [pqr_file, ion_strength, temperature, escl, vscl, threads, update_only_1],
                      [stage1_output, session_display])
        
        stage2_btn.click(run_stage2,
                      [n_reps, n_cycles, n_proteins, box_size, threads, sim_threads, update_only_2],
                      [stage2_output, session_display])
        
        stage3_btn.click(run_stage3,
                      [box_size, spacing, threads, update_only_3],
                      [stage3_output, session_display])
        
        pack_btn.click(pack_results, [skip_ve, skip_raw], [download_file, pack_status])
        
        plot_btn.click(create_plot,
                      [show_individual, show_pq, q_min, q_max, y_min, y_max, log_scale_x, log_scale_y],
                      [plot_image])

        # Update file downloads when plot is generated
        plot_btn.click(get_iqs_file, [], [iqs_file])
        plot_btn.click(get_stats_file, [], [stats_file])

        # Update jobs display when thread settings change
        threads.change(update_jobs_display, [threads, sim_threads], [jobs_display])
        sim_threads.change(update_jobs_display, [threads, sim_threads], [jobs_display])

        # Update concentration displays when parameters change
        n_proteins.change(update_concentrations, [n_proteins, box_size], 
                         [mw_display, conc_mg_ml_display, conc_mm_display])
        box_size.change(update_concentrations, [n_proteins, box_size], 
                       [mw_display, conc_mg_ml_display, conc_mm_display])

        # Update concentrations when Stage 1 completes
        stage1_btn.click(update_concentrations, [n_proteins, box_size], 
                      [mw_display, conc_mg_ml_display, conc_mm_display])

    return app

def main():
    """Launch the GUI"""
    app = create_gui()
    app.queue()
    app.launch(share=False, show_api=False, show_error=True)

if __name__ == "__main__":
    main()
