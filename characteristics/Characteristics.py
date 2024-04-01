import os
import pyBigWig
import math
import numpy as np
import torch
import psutil
from common import Chromosome_Info, Chromosome_Info_For_Simple
from sys import getsizeof
from characteristics.bwhelper import generate_array_from_matric, write_array_to_file, read_numbers_from_file

class Сharacteristic():

    def __init__(self, path :  str):
        self.path = path

class СharacteristicBigWig(Сharacteristic):

    def __init__(self, folder_res : str, size_of_bw, path : str = ""):
        super().__init__(
            path,
        )
        self.folder_res = folder_res
        self.size_of_bw = size_of_bw
        if (path != ""):
            self.prerpocess_all()
        print("before read bw")
        self.data_vector = []
        self.read_bw()
        self.chr_info = []
        with open(os.path.join(self.folder_res, "chr_info.pickle"), 'rb') as file:
            self.chr_info = pickle.load(file)
        
    
    def prerpocess_all(self):
        if not os.path.exists(self.folder_res):
            os.makedirs(self.folder_res)
        file_paths = []
        with open(self.path, 'r') as file:
            for line in file:
                file_paths.append(line.strip())      
        for file_name in file_paths:
            self.preprocess_bw(file_name)
    
    def preprocess_bw(self, file_name : str):
        print("current file: ", file_name)
        bw = pyBigWig.open(file_name)

        data_map = {}

        bw_size = bw.header()['nBasesCovered']
        chromosome_info_list = [] 
        j = 0
        for i in range(1, self.size_of_bw + 1):
            chr = "chr"
            if i == 23:
                chr += "X"
            elif i == 24:
                chr += "Y"
            else:
                chr += str(i)
            print(chr)

            chrom_size = bw.chroms(chr)

            values = bw.values(chr, 0, chrom_size)
            chromosome_info_list.append(Chromosome_Info(number=i, start_pos=j, lenght=chrom_size)) 

            for pos, value in enumerate(values, start=1):
                if not math.isnan(value):
                    data_map[j + pos] = value

            j += chrom_size
            del chrom_size
            del values
            del pos
        
        bw.close()

        # save hash table into file
        file_basename = os.path.basename(file_name)
        output_path = os.path.join(self.folder_res, file_basename.replace(".bw", "_m.pickle"))

        with open(output_path, "wb") as file:
            pickle.dump(data_map, file)

        chr_info_path = os.path.join(self.folder_res, "chr_info.pickle")
        # save chomosome info file
        with open(chr_info_path, 'wb') as file:
            pickle.dump(chromosome_info_list, file)


    def read_bw(self):
        print("read_bw")

        for filename in os.listdir(self.folder_res):
            if filename.endswith("m.pickle"):
                file_path = os.path.join(self.folder_res, filename)
                with open(file_path, "rb") as file:
                    loaded_map = pickle.load(file)
                    self.data_vector.append(loaded_map)

    def get_lines(self, chr : int, start : int, WINDOW_SIZE : int):
        result_tensor = torch.zeros((len(self.data_vector), WINDOW_SIZE))
        for idx, cur_map in enumerate(self.data_vector):
            new_start = self.chr_info[chr].start_pos + start
            for i in range(WINDOW_SIZE):
                if (new_start+i) in cur_map:
                    cur_map[new_start+i]
                    result_tensor[idx][i]
        return result_tensor

    def get_name(self):
        return "bw"

    def get_chr_info(self):
        return self.chr_info


class СharacteristicBigWigCSR(Сharacteristic):

    def __init__(self, folder_res : str, size_of_bw, path : str = ""):
        super().__init__(
            path,
        )
        self.folder_res = folder_res
        self.size_of_bw = size_of_bw
        self.index_in_row = []
        if (path != ""):
            self.prerpocess_all()
        if (len(self.index_in_row) == 0):
            self.index_in_row = np.load('resfold/my_array.npy')
        self.chr_info = [Chromosome_Info(number=0, start_pos=0, lenght=1000000)]

    def get_max_array_size(self):
        available_memory = psutil.virtual_memory().available
        dtype_size = np.dtype(np.int64).itemsize 
        max_elements = available_memory // dtype_size // self.big_wig_count
        # добавить деление на количество биг вигов
        return max_elements // 2
        
    
    def prerpocess_all(self):
        if not os.path.exists(self.folder_res):
            os.makedirs(self.folder_res)
        file_paths = []
        with open(self.path, 'r') as file:
            for line in file:
                file_paths.append(line.strip()) 
        self.big_wig_count = len(file_paths)   
        self.preprocess_subsequence(file_paths)   

    def preprocess_subsequence(self, file_paths):
        last_post_in_file = 0
        cur_pos = 0
        # я бы создала перед этим chr_info
        chrom_size = 1

        while cur_pos < chrom_size:
            cur_array = []
            cur_range = self.get_max_array_size() // 10
            for file_name in file_paths:
                print(file_name)
                bw = pyBigWig.open(file_name)
                chrom_size = bw.chroms('chr1')
                values = bw.values("chr1", cur_pos, min(cur_pos+cur_range, chrom_size))
                cur_array.append(values)
                del values
                del bw
            
            res_matrix = np.zeros((len(cur_array[0]), len(cur_array)), dtype=np.int32)
            for i in range(len(cur_array)):
                for j in range(len(cur_array[i])):
                    if not math.isnan(cur_array[i][j]):
                        res_matrix[j][i] = int(cur_array[i][j])
            del cur_array

            arr_col_data, index_in_row = generate_array_from_matric(res_matrix, last_post_in_file)
            del res_matrix
            last_post_in_file = index_in_row[len(index_in_row) - 1] // 2
            if cur_pos != 0:
                prev_index_row = np.load('resfold/my_array.npy')[:-1]
            else:
                prev_index_row = np.empty(0)
            index_in_row = np.append(prev_index_row, index_in_row)

            np.save('resfold/my_array.npy', index_in_row)
            self.index_in_row = index_in_row
            del index_in_row
            del prev_index_row
            if (cur_pos == 0):
                write_array_to_file(b'resfold/res.bin', np.array(arr_col_data, dtype=np.int64), b'wb')
            else:
                write_array_to_file(b'resfold/res.bin', np.array(arr_col_data, dtype=np.int64), b'a')
            cur_pos += cur_range

    
    def get_lines(self, chr : int, start : int, WINDOW_SIZE : int):
        start_in_file = self.index_in_row[start]
        end_posit = self.index_in_row[start + WINDOW_SIZE]
        answer = torch.zeros(10, WINDOW_SIZE)
        if (end_posit - start_in_file == 0):
            return answer

        cur_res = read_numbers_from_file('resfold/res.bin', end_posit - start_in_file, start_in_file)
        j = 0
        i = 0
        answer = torch.zeros(10, WINDOW_SIZE).type(torch.LongTensor)
        while j < len(self.index_in_row) and i < WINDOW_SIZE:
            count_in_cur_row = self.index_in_row[start + i + 1] - self.index_in_row[start + i]
            j_in_row = 0
            while (j_in_row < count_in_cur_row):
                answer[cur_res[j_in_row + j]][i] = cur_res[j_in_row + j + 1]
                j_in_row += 2
            j += j_in_row
            i += 1

        return answer

    def get_name(self):
        return "bwcsr"

    def get_chr_info(self):
        return self.chr_info


if __name__ == '__main__':
    chr_big_wig = СharacteristicBigWigCSR("/home/ojpochemy/dna-loader/resfold", 1, "/home/ojpochemy/dna-loader/file_paths.txt")
    # print(chr_big_wig.get_lines(0, 0, 10))
    # print(chr_big_wig.get_lines(0, 100, 10))
    print(chr_big_wig.get_name())

    



