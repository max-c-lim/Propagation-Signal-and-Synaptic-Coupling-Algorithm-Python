# Propagation-Signal-and-Synaptic-Coupling-Algorithm-Python

Python version of [Propagation-Signal-and-Synaptic-Coupling-Algorithm](https://github.com/ZhuoweiCheng/Propagation-Signal-and-Synaptic-Coupling-Algorithm)

REQUIRED PYTHON PACKAGES
numpy
multiprocessing

Notable changes made when converting from MATLAB to Python.

1. Instead of using [], use None instead.
2. Instead of 1 x N cell, use np.array with shape (N, )

CODE REQUIREMENTS/EXAMPLES
If on Windows, function must be with if __name__ == "__main__"

NOTES
In place operations could be messing things up

scan_reference_electrode results are different
rescan_candidate_cohorts results are different

Parallel processing can be optimized? See TODO

numpy operations can speed up calculations

Combine scan_reference_electrode and rescan_candidate_cohorts? (1 loop instead of 2)