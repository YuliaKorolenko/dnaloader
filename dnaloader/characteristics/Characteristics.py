import os
import pyBigWig
import math
import numpy as np
import torch
import psutil
import pickle
from dnaloader.common import Chromosome_Info, Chromosome_Info_For_Simple
from sys import getsizeof
from dataclasses import dataclass
from typing import List
from tqdm import tqdm
from memory_profiler import profile
from dnaloader.characteristics.bwhelper import write_array_to_file, read_numbers_from_file, read_and_convert_numbers_from_file, generate_array_from_tuples, create_upper_indixies, find_results, read_and_convert_light


@dataclass
class Bw_Meta:
    number_of_bw: int
    chromosome_info_list: List[Chromosome_Info]


class Characteristic():

    def __init__(self, path: str):
        self.path = path

    def get_chr_name(self, i: int):
        chr = "chr"
        if i == 23:
            chr += "X"
        elif i == 24:
            chr += "Y"
        else:
            chr += str(i)
        return (chr)


class CharacteristicBigWig(Characteristic):

    def __init__(self, folder_res: str, size_of_bw, path: str = ""):
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

    def preprocess_bw(self, file_name: str):
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
            chromosome_info_list.append(
                Chromosome_Info(
                    number=i,
                    start_pos=j,
                    lenght=chrom_size))

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
        output_path = os.path.join(
            self.folder_res, file_basename.replace(
                ".bw", "_m.pickle"))

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

    def get_lines(self, chr: int, start: int, WINDOW_SIZE: int):
        result_tensor = torch.zeros((len(self.data_vector), WINDOW_SIZE))
        for idx, cur_map in enumerate(self.data_vector):
            new_start = self.chr_info[chr].start_pos + start
            for i in range(WINDOW_SIZE):
                if (new_start + i) in cur_map:
                    cur_map[new_start + i]
                    result_tensor[idx][i]
        return result_tensor

    def get_name(self):
        return "bw"

    def get_chr_info(self):
        return self.chr_info


class CharacteristicBigWigCSR(Characteristic):

    # @profile
    def __init__(self, folder_res: str, size_of_bw: int,
                 type_of_loader: str = "hard", path: str = ""):
        super().__init__(
            path,
        )
        self.folder_res = folder_res
        self.size_of_bw = size_of_bw
        self.name_of_indexes = os.path.join(self.folder_res, 'indixies')
        self.type_of_loader = type_of_loader
        if (self.type_of_loader == "medium" or self.type_of_loader == "light"):
            self.name_of_indexes += ".npy"
        self.name_of_meta = os.path.join(self.folder_res, 'bw_meta.pickle')
        self.index_in_row = np.empty((0, 2), dtype=np.int64)
        if (path != ''):
            self.prerpocess_all()
        else:
            # проверить, что есть такие файлы, если нет, то кинуть ошибку и попросить перезагрузить
            # а тут в else флажок? и считывание из pybigWig
            # тут надо написать что-то типо, если хотите быстрее, то тынете
            # пропроцесс. Пока быстрых данных нет
            if (len(self.index_in_row) == 0 and self.type_of_loader == "medium" or
                    self.type_of_loader == "light"):
                self.index_in_row = np.load(self.name_of_indexes + '.npy')

            with open(self.name_of_meta, 'rb') as file:
                self.bw_meta = pickle.load(file)

    def get_max_array_size(self):
        available_memory = psutil.virtual_memory().available
        dtype_size = np.dtype(np.int64).itemsize
        max_elements = available_memory // dtype_size
        return max_elements // self.bw_meta.number_of_bw // 10

    # @profile
    def preprocess_meta(self, file_paths):
        chrom_lists = []
        file_name = file_paths[0]
        bw = pyBigWig.open(file_name)
        cur_pos = -1
        for i in range(1, self.size_of_bw + 1):
            cur_pos += 1
            chrom_size = bw.chroms(self.get_chr_name(i))
            chrom_lists.append(
                Chromosome_Info(
                    number=i,
                    start_pos=cur_pos,
                    lenght=chrom_size))
            cur_pos += chrom_size

        del bw, cur_pos, file_name, chrom_size
        # save metadata
        self.bw_meta = Bw_Meta(len(file_paths), chrom_lists)
        with open(self.name_of_meta, 'wb') as f:
            pickle.dump(self.bw_meta, f)

    def prerpocess_all(self):
        file_paths = []
        with open(self.path, 'r') as file:
            for line in file:
                file_paths.append(line.strip())
        # тут проверять, что есть хотя бы один биг виг
        if not os.path.exists(self.folder_res):
            os.makedirs(self.folder_res)
        self.preprocess_meta(file_paths)
        self.preprocess_subsequence(file_paths)

    def save_information(self, arr_col_data, chr_name: str,
                         cur_pos: int, index_in_row_cur, chr_num: int):
        mode = b'a'
        if (cur_pos == 0):
            mode = b'wb'
        if (self.type_of_loader == "medium"):
            self.index_in_row = np.append(
                self.index_in_row[:-1], index_in_row_cur)
            np.save(self.name_of_indexes + '.npy', self.index_in_row)
        elif (self.type_of_loader == "light"):
            self.__save__light_indexies(index_in_row_cur, cur_pos, chr_num)
        else:
            self.__save_hard_indexies(
                chr_name, index_in_row_cur, cur_pos, mode)
        del index_in_row_cur

        name_of_res_file = os.path.join(
            self.folder_res, "res_" + chr_name + ".bin")
        write_array_to_file(
            name_of_res_file,
            np.array(
                arr_col_data,
                dtype=np.int32),
            mode)

        del arr_col_data

    def __save_hard_indexies(
            self, chr_name: str, index_in_row_cur, cur_pos: int, mode):
        cur_name = self.name_of_indexes + "_" + chr_name + ".bin"
        if (cur_pos == 0):
            write_array_to_file(
                cur_name,
                np.array(
                    index_in_row_cur,
                    dtype=np.int32),
                mode)
        else:
            write_array_to_file(cur_name, np.array(
                index_in_row_cur[1:], dtype=np.int32), mode)

    def __save__light_indexies(self, index_in_row_cur, cur_pos: int, chr: int):
        # обычные индексы переводим в сжатые, прибавляем к существующим
        # результаты записываем в файл indexies
        upper_indexies = create_upper_indixies(
            np.array(index_in_row_cur), cur_pos + self.bw_meta.chromosome_info_list[chr - 1].start_pos)
        print("upper indexies: ", upper_indexies)
        self.index_in_row = np.concatenate(
            (self.index_in_row[:-1], upper_indexies), axis=0, dtype=np.int64)
        print(self.index_in_row)
        np.save(self.name_of_indexes + '.npy', self.index_in_row)

    def preprocess_subsequence(self, file_paths):
        for num_chr in range(1, self.size_of_bw + 1):
            last_post_in_file = 0
            cur_pos = 0
            chrom_size = self.bw_meta.chromosome_info_list[num_chr - 1].lenght
            chr_name = self.get_chr_name(num_chr)
            print("Chromosome preprocessing: ", chr_name)
            with tqdm(total=chrom_size, desc="Progress", unit="iter") as pbar:
                while cur_pos < chrom_size:
                    cur_array = []
                    cur_range = self.get_max_array_size()
                    for file_name in file_paths:
                        bw = pyBigWig.open(file_name)
                        values = bw.values(
                            chr_name, cur_pos, min(
                                cur_pos + cur_range, chrom_size))
                        cur_array.append(values)
                        del values
                        del bw

                    res_array = []
                    # count of nucleotids
                    for j in range(len(cur_array[0])):
                        # count of big wig
                        for i in range(len(cur_array)):
                            if not math.isnan(cur_array[i][j]):
                                res_array.append((i, j, int(cur_array[i][j])))
                    del i, j, cur_array

                    arr_col_data, index_in_row_cur = generate_array_from_tuples(
                        res_array, cur_range, last_post_in_file)
                    del res_array
                    last_post_in_file = index_in_row_cur[len(
                        index_in_row_cur) - 1] // 2

                    self.save_information(
                        arr_col_data, chr_name, cur_pos, index_in_row_cur, num_chr)
                    pbar.update(cur_range)
                    cur_pos += cur_range
            del chrom_size, chr_name

    def get_file_name(self, chr: int):
        return os.path.join(self.folder_res, "res_" +
                            self.get_chr_name(chr + 1) + ".bin")

    def get_indixies(self, WINDOW_SIZE: int, start: int, chr: int):
        if (self.type_of_loader == "medium"):
            start_pos_chr = self.bw_meta.chromosome_info_list[chr].start_pos
            start_in_array = start_pos_chr + start
            start_in_file = self.index_in_row[start_in_array]
            end_posit = self.index_in_row[start_in_array + WINDOW_SIZE]
            cur_indixies = self.index_in_row[start_in_array:
                                             start_in_array + WINDOW_SIZE + 1]
            return cur_indixies, start_in_file, end_posit
        elif (self.type_of_loader == "light"):
            # print((self.index_in_row[100]))
            gl_start = self.bw_meta.chromosome_info_list[chr].start_pos
            cur_indexies = find_results(
                self.index_in_row,
                gl_start + start,
                gl_start + start + WINDOW_SIZE)
            # print("cur indixies ", cur_indexies)
            return cur_indexies, cur_indexies[0][1], cur_indexies[len(
                cur_indexies) - 1][1]
        else:
            cur_name = self.name_of_indexes + "_" + \
                self.get_chr_name(chr + 1) + ".bin"
            cur_indixies = read_numbers_from_file(
                cur_name, WINDOW_SIZE + 1, start)
            end_posit = cur_indixies[WINDOW_SIZE]
            start_in_file = cur_indixies[0]
            return cur_indixies, start_in_file, end_posit

    def get_lines(self, chr: int, start: int, WINDOW_SIZE: int):
        cur_indixies, start_in_file, end_posit = self.get_indixies(
            WINDOW_SIZE, start, chr)
        # print(start_in_file, " ", end_posit)
        cur_res = read_numbers_from_file(
            self.get_file_name(chr),
            end_posit - start_in_file,
            start_in_file)
        # print("cur res: ", cur_res)
        if (self.type_of_loader == "light"):
            start_for_chr = self.bw_meta.chromosome_info_list[chr].start_pos
            # print(cur_res)
            # print(cur_indixies)
            return read_and_convert_light(WINDOW_SIZE, np.array(
                cur_res, dtype=np.int64), start_for_chr, cur_indixies, self.bw_meta.number_of_bw)
        return read_and_convert_numbers_from_file(WINDOW_SIZE, np.array(cur_res, dtype=np.int64), self.bw_meta.chromosome_info_list[chr].lenght,
                                                  np.array(cur_indixies, dtype=np.double), self.bw_meta.number_of_bw)

    def get_name(self):
        return "bwcsr"

    def get_chr_info(self):
        return self.bw_meta.chromosome_info_list


if __name__ == '__main__':
    chr_big_wig = CharacteristicBigWigCSR(
        "/home/ojpochemy/dna-loader/resfold",
        1,
        "/home/ojpochemy/dna-loader/file_paths.txt")
    # print(chr_big_wig.get_lines(0, 0, 10))
    # print(chr_big_wig.get_lines(0, 100, 10))
    print(chr_big_wig.get_name())
