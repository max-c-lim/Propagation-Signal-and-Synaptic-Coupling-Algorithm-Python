# Propagation-Signal-and-Synaptic-Coupling-Algorithm-Python

Python version of [Propagation-Signal-and-Synaptic-Coupling-Algorithm](https://github.com/ZhuoweiCheng/Propagation-Signal-and-Synaptic-Coupling-Algorithm)

## Installation
Download the file [automated_detection_propagation.py](automated_detection_propagation.py) and put it in your project folder, so that the function automated_detection_propagation can be imported.

## Required Python Libraries
1. numpy
2. pandas
3. multiprocessing

## Usage
The input of the algorithm is a np.array with shape (N,) where each element represents one electrode.
Each element is a np.array with shape (m,) that contains the spike times for each electrode.

Next, use the function automated_detection_propagation in the script automated_detection_propagation.py to extract the propagation signals in an array.
```python
from automated_detection_propagation import automated_detection_propagation

if __name__ == "__main__":
    list_of_propagation, propagating_times = automated_detection_propagation(spike_times, thres_freq, seconds_recording, thres_number_spikes, ratio, thres_cooccurrences, p)`
```
**Note**: The function `automated_detection_propagation` must be called somewhere within `if __name__ == "__main__":`

Use (`thres_freq` and `seconds_recording`) or `thres_number_spikes`:

If using (`thres_freq` and `seconds_recording`), use `None` for `thres_number_spikes`
```python
list_of_propagation, propagating_times = automated_detection_propagation(spike_times, thres_freq, seconds_recording, None, ratio, thres_cooccurrences, p)
```

If using `thres_number_spikes`, use `None` for `thres_freq` and `seconds_recording`
```python
list_of_propagation, propagating_times = automated_detection_propagation(spike_times, None, None, thres_number_spikes, ratio, thres_cooccurrences, p)
```

The parameters `small_window`, `big_window`, `ccg_before`, and `ccg_after` can be given input values in the function call
to alter how the algorithm generates the propagations. These four parameters and spike_times must be given in the same time unit
i.e. if spike_times contains times in ms, small_window must also be in ms. 
<br>**NOTE:** The default parameters are in ms, so if spike_times is not in ms, these values must be changed in the function call.
<br>The parameter `ccg_n_bins` can be changed as well.
```python
list_of_propagation, time_all = automated_detection_propagation(spike_times, None, None, 180, 0.5, 50, 50,
                                                                small_window=1.0, big_window=2.2,
                                                                ccg_before=2.0, ccg_after=2.0, ccg_n_bins=81)
```


Examples are provided in [example/example.py](example/example.py)

### Inputs to automated_detection_propagation
**spike_times:** np.array with shape (N, )
<br>
N columns represent N electrodes.
Each column contains a np.array with shape (m,) representing
the spike times for each electrode.

**thres_freq:** int, float, or None
<br>
A value representing the frequency lower bound of the spiking
frequency for all electrodes. Only electrodes that's above the
threshold will be considered as a reference electrode. For
example, enter 1 for 1Hz.

**seconds_recording:** int, float, or None
<br>
The length of recording in seconds. For example, 120 for
2 minutes recording.

**thres_number_spikes:** int, float, or None <br>
The threshold for total number of spikes on each reference
electrode. If is None, use thres_freq for threshold instead.

**ratio:** float
<br>
(Assuming the time unit is ms) let n1 denote the largest sum of counts in any 0.5 ms moving
    window in the cross-correlogram (CCG) and n2 denote the sum
    of counts of the 2 ms window with the location of the largest
    sum in the center. If the largest sum is found in the first
    1 ms or the last 1 ms of the CCG, take the sum of the counts
    of the first 2 ms window or the counts of the last 2 ms window
    as n2. This ratio is the lower bound threshold for n2/n1

**thres_cooccurrences:** int or float <br>
    Lower bound of the number of short latency co-occurrences each
    electrode needs to have.

**p:** int or float <br>
    Percentage of the maximum number of co-occurrences required for
    all constituent electrodes. p should be between 0 and 100.

#### Optional Inputs #### 
It is recommended to maintain the original values because these values
have been thoroughly tested to produce good results.

**small_window**: int or float: default=0.5 <br>
The size of the window of n1 described in ratio

**big_window**: int or float: default=2<br>
The size of the window of n2 described in ratio

**ccg_before** int or float: default=1.5<br>
The time before each reference spike to use when creating the CCG

**ccg_after** int or float: default=1.5<br>
The time after each reference spike to use when creating the CCG

**ccg_n_bins** int: default=61<br>
The number of bins to use when creating the CCG

### Outputs of automated_detection_propagation
**list_of_propagation:** np.array with shape (P,)<br>
    Each element contains a pandas.DataFrames (each with 4 columns) of electrode cohorts for
    each propagation (p) in a recording. Each DataFrame provides a list of candidate
    electrodes along with the latency between each electrode
    with the reference electrode, the number of co-occurrences,
    and the n2/n1 ratio.

**propagating_times:** np.array with shape (P,)<br>
Each element contains a 1d np.array with shape (Q,)
of spike times in the propagation with different number of anchor points chosen
for each propagation in list_of_propagation. The pth element in propagating_times
contains the spike times for the pth element in list_of_propagation.

The qth element in the inner array with shape (Q,) is an array containing the
propagating spike times isolated with q+1 anchor points. I.e. the 0th element
contains the propagating spike times isolated with 1 anchor point (an empty array),
the 1st element contains propagating spike times isolated with 2 anchor points,
the 2nd element contains propagating spike times isolated with 3 anchor points,
etc., until all constituent electrodes are used as anchor points.

The values in propagating_times will be in the same units as spike_times

## Notable changes made when converting from MATLAB to Python
1. Instead of using [], use None
2. Every instance of a 1 x N cell has been changed to a 1d np.array with shape (N, )
3. automated_detection_propagation is somewhat slower in Python than MATLAB
