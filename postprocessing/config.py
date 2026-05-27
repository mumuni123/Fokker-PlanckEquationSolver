from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

# Simulation output directory that contains density_*.dat and fields_*.dat.
DATA_DIR = ROOT / "output"

# All postprocessing figures are written here.
RESULTS_DIR = ROOT / "results"

# Source parameter file used to infer dt_snapshot, t_end, and electric-field normalization.
PARAMETERS_FILE = ROOT / "src" / "parameters.h"

# Settings for plot_time_evolution.py.
TIME_PREFIX = "fields"
POSITION_UM = 2.5
TIME_COLUMNS = ["Ex[V/m]"]

# Settings for plot_space_profile.py.
SPACE_FILE = DATA_DIR / "fields_00010.dat"
SPACE_COLUMNS = ["Ex[V/m]"]

# Plot style and output quality.
FIGSIZE = (7.2, 4.5)
DPI = 200
