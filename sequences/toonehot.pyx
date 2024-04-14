from libc.stdio cimport FILE, fopen, fclose, fread, fgets, fwrite, getline, fseek
from libc.stdlib cimport malloc, free, atoi, calloc
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

@cython.boundscheck(False)
@cython.wraparound(False) 
def read_n_values_from_dna(str filename, int window_size, long i, np.ndarray[np.int64_t, ndim=1] arr):
    cdef FILE *file
    cdef char *buffer = <char *>malloc(window_size)
    
    file = fopen(filename.encode(), "rb")
    if file == NULL:
        return "Error opening file"
    
    fseek(file, i, 0)
    fread(buffer, 1, window_size, file)
    
    fclose(file)

    cdef int sequence_len = len(buffer)
    cdef np.ndarray[FDTYPE_t, ndim=2] encoding = np.zeros((4, sequence_len), dtype=np.float32)
    cdef int j
    cdef char base

    for j in range(sequence_len):
        base = buffer[j]
        if base == ord('A'):
            encoding[arr[0], j] = 1
        elif base == ord('C'):
            encoding[arr[1], j] = 1
        elif base == ord('T'):
            encoding[arr[2], j] = 1
        elif base == ord('G'):
            encoding[arr[3], j] = 1
    free(buffer)
    return encoding

@cython.boundscheck(False)
@cython.wraparound(False) 
def toonehotbytes(bytes sequence, np.ndarray[np.int64_t, ndim=1] arr):
    cdef int sequence_len = len(sequence)
    cdef np.ndarray[FDTYPE_t, ndim=2] encoding = np.zeros((4, sequence_len), dtype=np.float32)
    cdef int j
    cdef char base

    for j in range(sequence_len):
        base = sequence[j]
        if base == ord('A'):
            encoding[arr[0], j] = 1
        elif base == ord('C'):
            encoding[arr[1], j] = 1
        elif base == ord('T'):
            encoding[arr[2], j] = 1
        elif base == ord('G'):
            encoding[arr[3], j] = 1
    return encoding