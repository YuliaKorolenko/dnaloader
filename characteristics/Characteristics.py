import os
import pyBigWig
import math
import numpy as np
import torch

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

            for pos, value in enumerate(values, start=1):
                if not math.isnan(value):
                    data_map[j + pos] = value

            j += chrom_size
        
        bw.close()

        # save hash table into file
        file_basename = os.path.basename(file_name)
        output_path = os.path.join(self.folder_res, file_basename.replace(".bw", "_hm.txt"))
        with open(output_path, "w") as file:
            for key, value in data_map.items():
                file.write(str(key) + "\t" + str(value) + "\n")

    def read_bw(self):
        print("read_bw")

        for filename in os.listdir(self.folder_res):
            data_map = {}
            
            file_path = os.path.join(self.folder_res, filename)
            with open(file_path, "r") as file:
                for line in file:
                    key, value = line.strip().split("\t")
                    data_map[key] = float(value)
            
            self.data_vector.append(data_map)

    def get_lines(self, chr : int, start : int, WINDOW_SIZE : int):
        result_tensor = torch.zeros((len(self.data_vector), WINDOW_SIZE))
        for idx, cur_map in enumerate(self.data_vector):
            # ans = torch.zeros(WINDOW_SIZE)
            # TODO: сделать не только для 1 хромосомы
            new_start = 0 + start
            for i in (0, WINDOW_SIZE):
                if (new_start+i) in cur_map:
                    result_tensor[idx][i] = cur_map[new_start+i]

        return result_tensor

    def get_name(self):
        return "bw"



if __name__ == '__main__':
    chr_big_wig = СharacteristicBigWig("/home/ojpochemy/dna-loader/resfold", 1)
    # print(chr_big_wig.get_lines(0, 0, 10))
    # print(chr_big_wig.get_lines(0, 100, 10))
    print(chr_big_wig.get_name())

    



