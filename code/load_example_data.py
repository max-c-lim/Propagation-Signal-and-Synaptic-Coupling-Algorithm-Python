# Load .mat file for comparison
example_data_path = "comparison_data/example_data.mat"

from scipy.io import loadmat
import numpy as np

THRES_FREQ = 1
SECONDS_RECORDING = 180
THRES_NUMBER_SPIKES = None  # 180
RATIO = 0.5
THRES_COOCCURRENCES = 50
P = 50

_example_data = loadmat(example_data_path)["spike_times"][0, :]
SPIKE_TIMES = np.array([a[0, :] if a.size > 0 else np.array([]) for a in _example_data], dtype=object)
