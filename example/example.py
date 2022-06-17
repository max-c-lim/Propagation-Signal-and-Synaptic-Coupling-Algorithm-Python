from automated_detection_propagation import automated_detection_propagation
from scipy.io import loadmat
import numpy as np

if __name__ == "__main__":
    """
    Load example_data.mat which contains example data from
    https://github.com/max-c-lim/Propagation-Signal-and-Synaptic-Coupling-Algorithm-Python/blob/main/example/example_data.mat
    """
    example_data_path = "example_data.mat"
    spike_times_raw = loadmat(example_data_path)["spike_times"][0, :]
    # Data must be reformatted
    spike_times = np.array([a[0, :] if a.size > 0 else np.array([]) for a in spike_times_raw], dtype=object)


    """
    spike_times is the spike times from 120 electrodes. Each element contains the
    spike time from one electrode. The time unit is ms. This recording is 180s in
    total, and 1Hz is used as the thresholding for the spiking frequency of the
    reference electrode. Put None for thres_number_spikes.
    """
    list_of_propagation, time_all = automated_detection_propagation(spike_times, 1, 180, None, 0.5, 50, 50)


    """
    Instead of thresholding for spiking frequency of the reference electrode,
    you can also use an absolute number to threshold the total number of
    spikes needed on each electrode. Put None for thres_freq and seconds_recording.
    Add the number you want to use for thres_number_spikes. Here, 180 is used.
    """
    list_of_propagation, time_all = automated_detection_propagation(spike_times, None, None, 180, 0.5, 50, 50)
