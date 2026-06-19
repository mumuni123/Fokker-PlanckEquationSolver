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
POSITION_UM = 1
TIME_COLUMNS = ["Ex[V/m]"]
TIME_X_AXIS_RANGE = None

# Settings for plot_space_profile.py.
SPACE_FILE = DATA_DIR / "fields_00020.dat"
SPACE_COLUMNS = ["Ex[V/m]"]
# Use None to plot the full available range, or set (min, max) to restrict it.
# Bounds can also be None individually, e.g. (1.0, None).
SPACE_X_AXIS_RANGE = None

# Settings for plot_background_density_evolution.py.
BACKGROUND_DENSITY_FILES = [
    DATA_DIR / "density_00035.dat",
]
BACKGROUND_DENSITY_COLUMN = "n_bkg_e[m^-3]"
BACKGROUND_DENSITY_X_AXIS_RANGE = None

# Settings for plot_beam_density_profile.py.
DENSITY_FILE = DATA_DIR / "density_00035.dat"
DENSITY_COLUMNS = ["n_beam[m^-3]"]
DENSITY_X_AXIS_RANGE = None

# Settings for plot_fv_profile.py.
FV_FILE = DATA_DIR / "fv_bkg_e_00005.dat"
FV_COLUMNS = ["F(u)"]
FV_X_AXIS_RANGE = (0.0, 3.0e8)

# Settings for plot_current_density_profile.py.
FACE_CURRENT_FILES = DATA_DIR / "current_00005.dat"
# Use None to plot every current column in FACE_CURRENT_FILES.
# Column names may omit units, e.g. "J_total_smoothed" matches "J_total_smoothed[A/m2]".
FACE_CURRENT_COLUMNS = [
    "J_total_face",
    #"J_bkg_e",
    #"J_beam",
]
FACE_CURRENT_X_AXIS_RANGE = None

# Plot style and output quality.
FIGSIZE = (7.2, 4.5)
DPI = 200
