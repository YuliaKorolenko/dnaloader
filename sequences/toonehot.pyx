import numpy as np

cimport cython
cimport numpy as np

ctypedef np.float32_t FDTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False) 
def toonehot(str sequence, dict base_to_index):
    cdef int sequence_len = len(sequence)
    cdef np.ndarray[FDTYPE_t, ndim=2] encoding = np.zeros((4, sequence_len), dtype=np.float32)
    cdef int j
    cdef str base

    for j in range(sequence_len):
        base = sequence[j]
        if base in base_to_index:
            encoding[base_to_index[base], j] = 1
    return encoding