cimport numpy as np 
import numpy as np 
from libc.stdlib cimport malloc, free 
from cpython cimport array 
from libc.stdio cimport FILE, fopen, fclose, fread, fgets, fwrite, getline, fseek
from libc.stdlib cimport malloc, free, atoi
cimport libc.stdlib
from array import array

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

def write_array_to_file(const char* filename, np.ndarray[np.int64_t, ndim=1] arr, mode):
    cdef FILE* file = fopen(filename, mode)  
    # Open file in binary write mode    
    if file == NULL:        
        raise FileNotFoundError("Failed to open file for writing")
    cdef int* ptr = <int*>libc.stdlib.calloc(len(arr), sizeof(int))
    for i in range(len(arr)):
        ptr[i] = arr[i]
    cdef size_t num_written = fwrite(ptr, sizeof(int), len(arr), file)
    fclose(file)
    return num_written

def read_numbers_from_file(filename, int n, int start_el):
    cdef char buffer[80000]
    cdef int number
    
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

    # идем по массиву, вычисляем сколько на каждую строку элементов, и потом их бахаем в answer из 
    # j - указатель в numbers
    # i - указатель в window_size
    j = 0

    fclose(file)
    return numbers