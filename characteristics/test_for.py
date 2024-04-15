from bwhelper import create_upper_indixies, find_results
import numpy as np


# aaa = create_upper_indixies(np.array([0, 0, 0, 3, 3, 4, 5, 5, 5, 6, 6, 6]), 0)
aaa = np.array([[0, 0],
                [3, 3],
                [5, 4],
                [6, 5],
                [9, 6]], dtype=np.int32)
index_in_row = np.empty((0, 2), dtype=np.int32)

new_array = np.concatenate((index_in_row, aaa[:-1]), axis=0, dtype=np.int32)
print(new_array)
