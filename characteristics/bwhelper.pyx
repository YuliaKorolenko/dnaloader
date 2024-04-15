cimport numpy as np 
import numpy as np 
from cpython cimport array 
from libc.stdio cimport FILE, fopen, fclose, fread, fgets, fwrite, getline, fseek
from libc.stdlib cimport malloc, free, atoi, calloc
cimport libc.stdlib
from array import array
import cython
from cython.parallel import prange

ctypedef np.int32_t FDTYPE_t

def generate_array_from_matric(np.ndarray[np.int32_t, ndim=2] arr, int start_pos):
    cdef int i, j, count_in_row
    
    arr_col_data = []
    index_pointers = []
    count_in_row = start_pos
    for j in range(arr.shape[0]):
        index_pointers.append(count_in_row * 2)
        for i in range(arr.shape[1]):
            if (arr[j, i] != 0):
                count_in_row += 1
                arr_col_data.append(i)
                arr_col_data.append(arr[j, i])
    index_pointers.append(count_in_row * 2)

    return arr_col_data, index_pointers

def generate_hic_from_tuples(list tuple_array, int nucleotids_size, int start_pos, int start_index):
    cdef int i, j, count_in_row, boarder
    arr_col_data = []
    index_pointers = []

    count_in_row = start_pos

    # row - (nucleotids number), coulumn - (difference), data
    j = 0
    for i in range(0, nucleotids_size):
        index_pointers.append(count_in_row * 2)
        boarder = start_index + i
        if (j < len(tuple_array)):
            while (j < len(tuple_array) and tuple_array[j][0] == boarder):
                count_in_row += 1
                arr_col_data.append(tuple_array[j][1] - tuple_array[j][0])
                arr_col_data.append(tuple_array[j][2])
                j += 1
    index_pointers.append(count_in_row * 2)
    del tuple_array
    return arr_col_data, index_pointers

def generate_array_from_tuples(list tuple_array, int nucleotids_size, int start_pos):
    cdef int i, j, count_in_row
    arr_col_data = []
    index_pointers = []

    count_in_row = start_pos

    # row - (big wig number), coulumn - (nucleotids number), data
    j = 0
    for i in range(0, nucleotids_size):
        index_pointers.append(count_in_row * 2)
        if (j < len(tuple_array)):
            while (j < len(tuple_array) and tuple_array[j][1] == i):
                count_in_row += 1
                arr_col_data.append(tuple_array[j][0])
                arr_col_data.append(tuple_array[j][2])
                j += 1
    index_pointers.append(count_in_row * 2)
    del tuple_array
    return arr_col_data, index_pointers

def write_array_to_file(filename, np.ndarray[np.int32_t, ndim=1] arr, mode):
    cdef FILE* file = fopen(filename.encode(), mode)  
    # Open file in binary write mode    
    if file == NULL:        
        raise FileNotFoundError("Failed to open file for writing")
    cdef int* ptr = <int*>libc.stdlib.calloc(len(arr), sizeof(int))
    for i in range(len(arr)):
        ptr[i] = arr[i]
    cdef size_t num_written = fwrite(ptr, sizeof(int), len(arr), file)
    fclose(file)
    return num_written

def read_numbers_from_file(filename, long n, long start_el):
    #cdef char buffer[10000000]
    cdef unsigned char* buffer = <unsigned char*>(calloc(n, sizeof(int)));
    cdef unsigned int number
    
    cdef FILE* file = fopen(filename.encode(), "rb")
    if file == NULL:
        print("Error opening file")
        return None
    fseek(file, start_el * sizeof(int), 0)
    numbers = []

    cdef int j
    fread(buffer, sizeof(int), n, file)
    for i in range(0, n):
        j = i * 4
        number = buffer[j + 0] + (buffer[j + 1] << 8) + (buffer[j + 2] << 16) + (buffer[j + 3] << 24)
        numbers.append(number)
    free(buffer)

    # идем по массиву, вычисляем сколько на каждую строку элементов
    # j - указатель в numbers
    # i - указатель в window_size

    fclose(file)
    return numbers

@cython.boundscheck(False)
@cython.wraparound(False)
def read_and_convert_numbers_from_file(int window_size, np.ndarray[np.int64_t, ndim=1] cur_res, int len_chr,
                                       np.ndarray[np.double_t, ndim=1] indixies, int bw_count):
    cdef np.ndarray[FDTYPE_t, ndim=2] encoding = np.zeros((bw_count, window_size), dtype=np.int32)
    if (len(cur_res) == 0):
        return encoding
    cdef int j = 0
    cdef int i = 0
    cdef int count_in_cur_row, j_in_row
    while j < len_chr and i < window_size:
        count_in_cur_row = int(indixies[i + 1] - indixies[i])
        j_in_row = 0
        while (j_in_row < count_in_cur_row):
            encoding[cur_res[j_in_row + j]][i] = cur_res[j_in_row + j + 1]
            j_in_row += 2
        j += j_in_row
        i += 1
    return encoding

@cython.boundscheck(False)
@cython.wraparound(False)
def read_and_convert_light(int window_size, np.ndarray[np.int64_t, ndim=1] cur_res, long len_chr,
                                       np.ndarray[np.int64_t, ndim=2] up_indx, int bw_count):
    cdef np.ndarray[FDTYPE_t, ndim=2] encoding = np.zeros((bw_count, window_size), dtype=np.int32)
    if (len(cur_res) == 0):
        return encoding
    cdef int j = 0
    cdef int i = 0
    cdef int cur_up = 0
    cdef long start = up_indx[0][0]
    cdef int count_in_cur_row, j_in_row
    while j < len_chr and i < window_size:
        count_in_cur_row = 0
        if (cur_up + 1 < len(up_indx) and i + start == up_indx[cur_up + 1][0]):
            count_in_cur_row = int(up_indx[cur_up + 1][1] - up_indx[cur_up][1])
            cur_up += 1
        j_in_row = 0
        while (j_in_row < count_in_cur_row):
            encoding[cur_res[j_in_row + j]][i] = cur_res[j_in_row + j + 1]
            j_in_row += 2
        j += j_in_row
        i += 1
    return encoding

@cython.boundscheck(False)
@cython.wraparound(False)
def read_convert_to_hic_matrix(str file_indexies, str file_vals, int start, int window_size):
    indixies = read_numbers_from_file(file_indexies, window_size + 1, start)
    end_pos = indixies[window_size]
    start_in_file = indixies[0]
    
    cur_res = read_numbers_from_file(file_vals, end_pos - start_in_file, start_in_file)

    cdef np.ndarray[FDTYPE_t, ndim=2] hic_matrix = np.zeros((window_size, window_size), dtype=np.int32)
    if (len(cur_res) == 0):
        return hic_matrix
    cdef int ofset
    cdef int j = 0
    cdef int i = 0
    cdef int count_in_cur_row, j_in_row
    while i < window_size:
        count_in_cur_row = int(indixies[i + 1] - indixies[i])
        j_in_row = 0
        while (j_in_row < count_in_cur_row):
            # j - отступ от текущей диагонали
            offset = i + cur_res[j_in_row + j]
            if (offset < window_size):
                hic_matrix[i][offset] = cur_res[j_in_row + j + 1]
                hic_matrix[offset][i] = cur_res[j_in_row + j + 1]
            j_in_row += 2
        j += j_in_row
        i += 1
    return hic_matrix

@cython.boundscheck(False)
@cython.wraparound(False)
def create_upper_indixies(np.ndarray[np.int64_t, ndim=1] indx, long cur_indx):
    cdef list result = []
    cdef int i = 0
    cdef long start

    while i < indx.shape[0]:
        start = i
        while i < indx.shape[0] - 1 and indx[i] == indx[i + 1]:
            i += 1
        if (cur_indx + start < 0):
            print(cur_indx, " ", start)
        result.append((cur_indx + start, indx[start]))
        i += 1

    return np.array(result, dtype=np.int64)

@cython.boundscheck(False)
@cython.wraparound(False)
def find_results(np.ndarray[np.int64_t, ndim=2] upper_indx, long start_val, long end_val):
    cdef int start_pos = binary_search(upper_indx, start_val)
    cdef int end_pos = binary_search(upper_indx, end_val)
    return upper_indx[start_pos : end_pos + 1]

@cython.boundscheck(False)
@cython.wraparound(False)
cdef binary_search(np.ndarray[np.int64_t, ndim=2] upper_indx, long k): 
    low = 0 
    high = len(upper_indx) - 1 
    result = -1
    while low <= high:
        mid = (low + high) // 2
        if upper_indx[mid][0] <= k:
            result = mid
            low = mid + 1
        else:
            high = mid - 1

    return result