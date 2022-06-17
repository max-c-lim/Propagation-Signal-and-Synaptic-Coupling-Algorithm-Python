import numpy as np
from pandas import DataFrame
from multiprocessing import Pool


def int_round(n):
    """
    Rounds a float (n) to the nearest integer
    Uses standard math rounding, i.e. 0.5 rounds up to 1
    """

    n_int = int(n)
    if n - n_int >= 0.5:
        n_int += 1
    return n_int


def get_propagation_time(list_of_propagation, spike_times):
    """
    This function generates a sequence of eAPs for each
    propagation using different number of anchor points

    Inputs:
        list_of_propagation: np.array
            Output of get_propagation with shape (P,)
        spike_times: np.array
            With shape (N,) where N elements represent N electrodes.
            Each column contains a np.array with shape (m,) representing
            the spike times for each electrode. The spike times should be
            in units of ms. This code deals with data sampled at 20000Hz.
            Upsample or downsample your data into 20000Hz sample rate before
            feeding into the algorithm.

    Output:
        propagating_times: np.array
            With shape (P,) where each element contains a 1d np.array of spike
            times in the propagation with different number of anchor points chosen
            for each propagation in list_of_propagation with the same order.
            The 0th element in each np.array of propagating_times is np.array
            containing the propagating spike times isolated with 2 anchor points,
            the 1st element is the propagating spike times isolated with 3 anchor points,
            etc., until all constituent electrodes are used as anchor points.
    """

    n = list_of_propagation.size
    propagating_times = []
    for i in range(n):
        time_signal = []
        time = []
        current_signal = list_of_propagation[i]

        ind0 = 0
        index = np.flatnonzero(current_signal.loc[:, "latency"] > 0.001)
        index = np.concatenate(([ind0], index))
        current_signal = current_signal.loc[index, :].copy()
        size_of_signal, _ = current_signal.shape

        current_signal.loc[ind0, "small_window_cooccurrences"] = -888
        ID0 = current_signal.loc[ind0, "ID"]
        reference_spike_times = spike_times[ID0]

        for n_ele in range(1, size_of_signal):  # ele = electrode
            ind1 = current_signal.loc[:, "small_window_cooccurrences"].idxmax()
            current_signal.loc[ind1, "small_window_cooccurrences"] = -888
            ID1 = current_signal.loc[ind1, "ID"]
            target_spike_times = spike_times[ID1]

            ind_all = []
            for j in range(reference_spike_times.size):
                ref = reference_spike_times[j]
                index = np.flatnonzero((target_spike_times > ref) & (target_spike_times < ref+1.5))
                if index.size > 0:
                    ind_all.append(j)

            time.extend(reference_spike_times[ind_all])
            time_signal.append(np.sort(np.unique(time)))
        propagating_times.append(np.asarray(time_signal, dtype=object))

    return np.asarray(propagating_times, dtype=object)


def get_propagation(electrode_cohorts):
    """
    This function generates a collection of cohort electrodes, each
    representing an eAP propagation in each recording

    Inputs:
        electrode_cohorts: np.array
            Output of rescan_candidate_cohorts

    Output:
        list_of_propagation: np.array
            Contains a pandas.DataFrame of electrode cohorts for each propagation
            in a recording. Each DataFrame provides a list of candidate
            electrodes along with the latency between each electrode
            with the reference electrode, the number of co-occurrences,
            and the n2/n1 ratio.
    """

    list_of_propagation = []
    for i in range(electrode_cohorts.size):
        if electrode_cohorts[i].size > 0:
            temps = electrode_cohorts[i]
            m1, n1 = temps.shape
            if temps.size > 0 and not np.any(temps[1, :] < 0) and n1 > 1:
                index = np.flatnonzero(temps[0, :] == i)
                if index.size > 0:
                    order = np.argsort(temps[1, :])
                    order = np.delete(order, np.flatnonzero(order == index))
                    order = np.concatenate((index, order))
                    temps = temps[:, order]

                    table_data = {
                        "ID": temps[0, :],
                        "latency": temps[1, :],
                        "small_window_cooccurrences": temps[2, :],
                        "n1_n2_ratio": temps[3, :],
                    }
                    table = DataFrame(table_data).astype({"ID": int})
                    if table.size > 0:
                        list_of_propagation.append(table)

    return np.asarray(list_of_propagation, dtype=object)


def rescan_candidate_cohorts(candidate_cohorts, thres_cooccurrences, p):
    """
    This function rescans each set of candidate electrodes found for each
    reference electrode. First, find the electrode with the maximum number of
    co-occurrences with the reference electrode. Second, scan through all
    other electrodes in the set of candidate electrodes, to identify only
    the electrodes with more than p * the maximum number of co-occurrences
    and more than thres_number_spikes co-occurrences in the 0.5ms window with
    maximum number of co-occurrences in the CCG. The electrodes that satisfy
    this criterion are kept in electrode_cohorts

    Inputs:
        candidate_cohorts: np.array
            Output of scan_reference_electrode.py. np.array contains a
            list of candidate constituent electrodes for each reference
            electrode. Each row provides a list of candidate electrodes
            along with the latency between each electrode with the
            reference electrode, the number of co-occurrences and the n2/n1
            ratio.
        thres_cooccurrences: int or float
            Lower bound of the number of short latency co-occurrences each
            electrode needs to have.
        p: int or float
            Percentage of the maximum number of co-occurrences required for
            all constituent electrodes. p should be between 0 and 100.

    Output:
        electrode_cohorts: np.array
            Contains the constituent electrodes for each reference electrode.
            Each element is a np.array that contains the candidate electrodes
            along with the latency between each electrode with the reference electrode,
            the number of cooccurrences, and the n2/n1 ratio.

    """

    numbers = candidate_cohorts.size
    electrode_cohorts = [np.array([])] * numbers

    for i in range(numbers):
        if candidate_cohorts[i].size > 0:
            current_cohort = candidate_cohorts[i]

            reference = np.flatnonzero(current_cohort[0, :] == i)

            target_electrodes = np.flatnonzero((current_cohort[2, :] >= thres_cooccurrences) & (current_cohort[1, :] != 0) & (current_cohort[0, :] != reference))

            if target_electrodes.size > 0:
                cooccurrences = p/100 * max(current_cohort[2, target_electrodes])
                index = np.flatnonzero(current_cohort[2, :] >= cooccurrences)
                index = np.union1d(index, reference)
                current_cohort_new = current_cohort[:, index]
                electrode_cohorts[i] = current_cohort_new

    return np.array(electrode_cohorts, dtype=object)


def scan_reference_electrode(spike_times, thres_freq, seconds_recording, thres_number_spikes, ratio):
    """
    This function generates a cell array containing candidate electrode cohorts
    for each electrode. Each cell corresponds to an electrode with the same
    order of your input. If a cell is empty, there's no electrode cohorts
    associated with this electrode. Use rescan_each_reference to find
    constituent electrodes from candidate electrode cohorts.

    Inputs:
        spike_times: np.array
            With shape (N,) where N columns represent N electrodes.
            Each column contains a np.array with shape (m,) representing
            the spike times for each electrode. The spike times should be
            in units of ms. This code deals with data sampled at 20000Hz.
            Upsample or downsample your data into 20000Hz sample rate before
            feeding into the algorithm.
        thres_freq: int or float
            A value representing the frequency lower bound of the spiking
            frequency for all electrodes. Only electrodes that's above the
            threshold will be considered as a reference electrode. For
            example, enter 1 for 1Hz.
        seconds_recording: int or float
            The length of recording in seconds. For example, 120 for
            2 minutes recording.
        thres_number_spikes: int or None
            The threshold for total number of spikes on each reference
            electrode. If is None, use thres_freq for threshold instead.
        ratio: float
            Let n1 denote the largest sum of counts in any 0.5 ms moving
            window in the cross-correlogram (CCG) and n2 denote the sum
            of counts of the 2 ms window with the location of the largest
            sum in the center. If the largest sum is found in the first
            1 ms or the last 1 ms of the CCG, take the sum of the counts
            of the first 2 ms window or the counts of the last 2 ms window
            as n2. This ratio is the lower bound threshold for n2/n1
    Output:
        candidate_cohorts: np.array
            Contains a list of candidate constituent electrodes for each reference
            electrode. Each row provides a list of candidate electrodes
            along with the latency between each electrode with the
            reference electrode, the number of co-occurrences and the n2/n1
            ratio.

    Notes:
        61 is for 61 bins in CCG
    """

    if thres_number_spikes is not None:
        thres = thres_number_spikes
    else:
        thres = thres_freq * seconds_recording

    n_e = spike_times.size
    small_window = 0.5
    t = int_round(small_window/0.05)
    candidate_cohorts = [np.array([])] * n_e

    # Used to replicate MATLAB's parfor's multiprocessing
    with Pool(initializer=_scan_reference_electrode_worker_init, initargs=(spike_times, thres, ratio, n_e, t)) as pool:
        for i, cohort_data in enumerate(pool.map(_scan_reference_electrode_func, range(n_e))):
            if cohort_data is not None:
                candidate_cohorts[i] = cohort_data

    return np.array(candidate_cohorts, dtype=object)


def _scan_reference_electrode_worker_init(spike_times, thres, ratio, n_e, t):
    # Initialize variables for parallel processing worker
    # TODO: COULD BE CONSUMING TOO MUCH RAM SINCE "spike_times" IS VERY LARGE AND COPIED FOR EACH WORKER. INVESTIGATE FURTHER. cache in memory?
    global _scan_reference_electrode_parfor_dict
    _scan_reference_electrode_parfor_dict = {
        "spike_times": spike_times,
        "thres": thres,
        "ratio": ratio,
        "n_e": n_e,
        "t": t,
    }


def _scan_reference_electrode_func(electrode):
    # Function that each parallel processing worker will execute

    spike_times = _scan_reference_electrode_parfor_dict["spike_times"]
    thres = _scan_reference_electrode_parfor_dict["thres"]
    ratio = _scan_reference_electrode_parfor_dict["ratio"]
    n_e = _scan_reference_electrode_parfor_dict["n_e"]
    t = _scan_reference_electrode_parfor_dict["t"]

    ref = spike_times[electrode]
    n = ref.size  # Number of spikes detected by electrode
    # Assume in ms scale
    if n >= 0 and n >= thres:
        time_delay = np.full(n_e, -999.0)
        small_window_cooccurrences = np.zeros(n_e)
        spikes_ratio = np.full(n_e, -999.0)
        for electrode2 in range(n_e):
            if electrode2 == electrode:
                continue

            ccg = np.zeros(61)  # ccg = crosscorrelogram
            tar = spike_times[electrode2]  # tar = target(s)
            for k in range(n):
                ref_value = ref[k]
                index = np.flatnonzero((tar >= ref_value-1.5) * (tar <= ref_value + 1.5))
                if index.size > 0:
                    for i in index:
                        bin_value = tar[i] - ref_value
                        ccg[int_round(bin_value / 0.05 + 30)] += 1  # + 30 since index 30 is the center of the 61-bin ccg
            spikes_small_window = ccg[:t+1].sum()
            location = 0
            for i in range(61-t):
                sum_small_window = ccg[i:i+t+1].sum()
                if sum_small_window > spikes_small_window:
                    spikes_small_window = sum_small_window
                    location = i

            delay = location + np.argmax(ccg[location:location+t+1])
            if 20 < delay < 40:
                spikes_big_window = ccg[delay - 20:delay + 20 + 1].sum()
            elif delay <= 20:
                spikes_big_window = ccg[:41].sum()
            else:  # delay >= 40:
                spikes_big_window = ccg[20:61].sum()

            if spikes_small_window >= ratio * spikes_big_window and spikes_big_window >= 1:
                time_delay[electrode2] = (delay - 30) * 0.05
                small_window_cooccurrences[electrode2] = spikes_small_window
                spikes_ratio[electrode2] = spikes_small_window / spikes_big_window

        # Set values for own electrode
        time_delay[electrode] = 0
        small_window_cooccurrences[electrode] = n
        spikes_ratio[electrode] = 1

        ind = np.flatnonzero(time_delay >= -1.5)  # ind = indices which are the electrode ids that are a part of the cohort
        return np.vstack((ind, time_delay[ind], small_window_cooccurrences[ind], spikes_ratio[ind]))


def automated_detection_propagation(spike_times, thres_freq, seconds_recording, thres_number_spikes, ratio, thres_cooccurrences, p):
    """
    This function detects eAP propagation in a recording and generates a
    collection of cohort electrodes, each representing an eAP propagation in
    each recording. Along with the cohort electrodes, this function also
    outputs the propagating spike times for each eAP propagation with
    different number of anchor points.

    Inputs:
        spike_times: np.array
            With shape (N,) where N columns represent N electrodes.
            Each column contains a np.array with shape (m,) representing
            the spike times for each electrode. The spike times should be
            in units of ms. This code deals with data sampled at 20000Hz.
            Upsample or downsample your data into 20000Hz sample rate before
            feeding into the algorithm.
        thres_freq: int, float, or None
            A value representing the frequency lower bound of the spiking
            frequency for all electrodes. Only electrodes that's above the
            threshold will be considered as a reference electrode. For
            example, enter 1 for 1Hz.
        seconds_recording: int, float, or None
            The length of recording in seconds. For example, 120 for
            2 minutes recording.
        thres_number_spikes: int, float, or None
            The threshold for total number of spikes on each reference
            electrode. If is None, use thres_freq for threshold instead.
        ratio: float
            Let n1 denote the largest sum of counts in any 0.5 ms moving
            window in the cross-correlogram (CCG) and n2 denote the sum
            of counts of the 2 ms window with the location of the largest
            sum in the center. If the largest sum is found in the first
            1 ms or the last 1 ms of the CCG, take the sum of the counts
            of the first 2 ms window or the counts of the last 2 ms window
            as n2. This ratio is the lower bound threshold for n2/n1
        thres_cooccurrences: int or float
            Lower bound of the number of short latency co-occurrences each
            electrode needs to have.
        p: int or float
            Percentage of the maximum number of co-occurrences required for
            all constituent electrodes. p should be between 0 and 100.

    Outputs:
        list_of_propagation: np.array
            Contains pandas.DataFrames of electrode cohorts for each propagation
            in a recording. Each DataFrame provides a list of candidate
            electrodes along with the latency between each electrode
            with the reference electrode, the number of co-occurrences,
            and the n2/n1 ratio.
        propagating_times: np.array
            With shape (P,) where each element contains a 1d np.array of spike
            times in the propagation with different number of anchor points chosen
            for each propagation in list_of_propagation with the same order.
            The 0th element in each np.array of propagating_times is np.array
            containing the propagating spike times isolated with 2 anchor points,
            the 1st element is the propagating spike times isolated with 3 anchor points,
            etc., until all constituent electrodes are used as anchor points.

    Possible optimizations:
        Optimize parallel processing. See TODO
        numpy operations can speed up calculations
        Combine scan_reference_electrode and rescan_candidate_cohorts? (1 loop instead of 2)
    """

    candidate_cohorts = scan_reference_electrode(spike_times, thres_freq, seconds_recording, thres_number_spikes, ratio)
    electrode_cohorts = rescan_candidate_cohorts(candidate_cohorts, thres_cooccurrences, p)
    list_of_propagation = get_propagation(electrode_cohorts)
    propagating_times = get_propagation_time(list_of_propagation, spike_times)
    return list_of_propagation, propagating_times
