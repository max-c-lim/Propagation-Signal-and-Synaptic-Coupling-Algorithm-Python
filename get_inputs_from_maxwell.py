import h5py
from collections import OrderedDict
import numpy as np
from scipy.io import savemat


def get_inputs_from_maxwell(path, rec_id="0000", well_id="000", first_n_min=None, only_routed_electrodes=True, raw_data_recorded=True):
    """
    Script to retrieve input data (spike_times) to automated_detection_propagation
    from raw Maxwell recording (.h5 format)

    Inputs:
        path: str
            Path to raw recording file
        rec_id: str
            Which recording to take data from if the recording file has multiple recordings
        well_id: str
            Which well to take data from if the recording file has multiple wells
        first_n_min: None or int
            Only the spike times from the first first_n_min of recording will be outputted
            If None, entire recording will be used
        only_routed_electrodes: bool
            If True, only routed electrodes' spike times will be returned
            If False, all electrodes' spike times will be returned
        raw_data_recorded: bool
            If True, time point 0 for the spike times will be set to the starting time of the recording. Any spike that occurs before the recording start time will be discarded
            If False, time point 0 will be the first spike time - 1

    Output:
        spike_times: list
            Input to automated_detection_propagation (in ms). Only electrodes with at least 1 spike will be returned
        channel_map: np.array
            1d array with same length as spike_times. Maps index in spike_times to channel id

    Notes:
        This function assumes that the spike times stored in the raw .h5 maxwell recording file are sorted
    """

    """
    Exploring structure of Maxwell recordings
    rec.visit(lambda x: print(x)) 
        assay
        assay/inputs
        assay/inputs/electrodes   # electrodes: list = json.loads(rec["assay/inputs/electrodes"][0])["electrodes"]["0"]
        assay/inputs/record_time  # record_time: int = int(rec["assay/inputs/record_time"][0])
        assay/inputs/spike_only
        assay/run_id
        assay/script_id
        bits
        data_store
        data_store/data0000
        data_store/data0000/events
        data_store/data0000/groups
        data_store/data0000/groups/routed   # Data in this path will not exist if recording is not recorded first 
        data_store/data0000/groups/routed/channels
        data_store/data0000/groups/routed/frame_nos  # Frame number of each sample in raw
        data_store/data0000/groups/routed/raw
        data_store/data0000/groups/routed/triggered
        data_store/data0000/recording_id
        data_store/data0000/settings
        data_store/data0000/settings/gain
        data_store/data0000/settings/hpf
        data_store/data0000/settings/lsb
        data_store/data0000/settings/mapping  # an array where each element is (channel, electrode, x, y)
        data_store/data0000/settings/sampling  # in Hz
        data_store/data0000/settings/spike_threshold
        data_store/data0000/spikes      # in samples (aka frames)
        data_store/data0000/start_time  # in milliseconds 
        data_store/data0000/stop_time   # in milliseconds
        data_store/data0000/well_id
        environment
        environment/diagnosis
        environment/temperature
        hdf_version
        mxw_version
        notes
        recordings
        recordings/rec0000  # /well000 contains HardLinks to data in data_store/data0000
        version
        wellplate
        wellplate/id
        wellplate/version
        wellplate/well000
        wellplate/well000/control
        wellplate/well000/group_color
        wellplate/well000/group_name
        wellplate/well000/id
        wellplate/well000/name
        wells
        wells/well000
        
    All data is stored in array (even singles numbers)
    """
    h5 = h5py.File(path)
    rec_key = f"recordings/rec{rec_id}/well{well_id}"
    rec = h5[rec_key]

    sampling_hz = rec["settings/sampling"][0]
    sampling_khz = sampling_hz / 1000
    spikes = rec["spikes"]

    if raw_data_recorded:
        frame_nos = rec["groups/routed/frame_nos"]
        first_frame = frame_nos[0]
        last_frame = frame_nos[-1]
    else:
        first_frame = spikes[0][0]
        last_frame = spikes[-1][0]

    if first_n_min is not None:
        first_n_samples = first_n_min * 60 * sampling_hz
        last_frame = min(last_frame, first_frame + first_n_samples)

    if only_routed_electrodes:
        channels_selected = set(rec["groups/routed/channels"])
    else:
        channels_selected = None
    channel_spikes = OrderedDict()
    channel_map = []

    for spike in spikes:
        frame_num = spike[0]
        if frame_num < first_frame:
            continue
        elif frame_num > last_frame:
            break

        spike_ms = (frame_num - first_frame) / sampling_khz
        spike_channel = spike[1]
        if channels_selected is None or spike_channel in channels_selected:
            if spike_channel not in channel_spikes:
                channel_spikes[spike_channel] = [spike_ms]
                channel_map.append(spike_channel)
            else:
                channel_spikes[spike_channel].append(spike_ms)

    spike_times = [np.asarray(st) for st in channel_spikes.values()]
    return spike_times, channel_map


if __name__ == "__main__":
    recording = "220707_16460_000439_data"
    spike_times, channel_map = get_inputs_from_maxwell(recording+".raw.h5", first_n_min=5)
    print(f"N spikes: {sum([len(spikes) for spikes in spike_times])}")
    spike_times_np = np.asarray(spike_times, dtype=object)
    np.save(f"{recording}_spike_times.npy", spike_times_np)
    savemat(f"code/matlab/{recording}_spike_times.mat", {"spike_times": spike_times_np})

