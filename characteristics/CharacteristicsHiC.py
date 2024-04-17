from .Characteristics import Characteristic
import cooler
import time
import numpy as np
import math
from characteristics.bwhelper import generate_array_from_tuples, write_array_to_file, read_numbers_from_file, read_convert_to_hic_matrix, generate_hic_from_tuples
import os
from dataclasses import dataclass
import pickle
import pandas
from typing import List
from common import Chromosome_Info
import h5py


@dataclass
class HiC_Meta:
    bin_size: int
    chromsizes : List[Chromosome_Info]


class CharacteristicHiCColer(Characteristic):
    def __init__(self, hic_path: str):
        super().__init__(
            hic_path,
        )
        self.c_matrix = cooler.Cooler(f"{hic_path}::/resolutions/1000")

    def get_lines(self, chr: int, start: int, WINDOW_SIZE: int):
        return self.c_matrix.matrix(balance=False).fetch(
            (self.get_chr_name(chr + 1), start, start + WINDOW_SIZE))

    def get_name(self):
        return "hi_c"


class CharacteristicHiCWithLimit(Characteristic):
    def __init__(self, hic_path: str, folder_res: str):
        super().__init__(
            hic_path,
        )
        self.folder_res = folder_res
        self.indexes_name = os.path.join(self.folder_res, "indexies.bin")
        self.vals_name = os.path.join(self.folder_res, "vals")
        self.meta_name = os.path.join(self.folder_res, "hic_meta.pickle")
        self.chrom_count = 22
        with open(self.meta_name, 'rb') as file:
            self.hic_meta = pickle.load(file)
            print("size meta: ", len(self.hic_meta.chromsizes))
        self.up_boader = 1_000_000

    def __preprocess_meta(self, chromsizes, bin_size : int) -> None:
        chrom_lists = []
        for i in range(0, self.chrom_count):
            # + chr так как в масиве для каждой хромосомы - indexies на 1 больше элемент, чем размер
            chrom_lists.append(
                Chromosome_Info(
                    number=i,
                    start_pos=chromsizes[i] + i,
                    lenght=0))

        print("chrom lists: ", chrom_lists)
        
        self.hic_meta = HiC_Meta(bin_size=bin_size, chromsizes=chrom_lists)
        with open(self.meta_name, 'wb') as file:
            pickle.dump(self.hic_meta, file)

    def __save_values(self, j : int, i : int, indexies, vals, chr_name : str) -> None:
        mode = b'a'
        if (j == 0):
            mode = b'wb'
        
        mode_for_index = mode
        if (i > 1):
            mode_for_index = b'a'

        if (j == 0):
            write_array_to_file(
                self.indexes_name, np.array(
                    indexies, dtype=np.int32), mode_for_index)
        else:
            write_array_to_file(self.indexes_name, np.array(
                indexies[1:], dtype=np.int32), mode_for_index)

        write_array_to_file(
            self.vals_name + chr_name + ".bin", np.array(
                vals, dtype=np.int32), mode)
        del indexies
        del vals

    def preprocess(self) -> None:
        if not os.path.exists(self.folder_res):
            os.mkdir(self.folder_res)
        cmat = cooler.Cooler(f"{self.path}::/resolutions/1000")

        hdf5_file = h5py.File(self.path, 'r')
        self.__preprocess_meta(hdf5_file['resolutions'][self.hic_meta.bin_size]['indexes']['chrom_offset'], cmat.info['bin-size'])

        pixels = cmat.pixels()
        start_time = time.time()
        batch_size = 30_000
        cur_bin = 0

        for i in range(1, self.chrom_count + 1):
            chr_name = self.get_chr_name(i)
            last_size = 0
            print(chr_name)
            bin_for_chr = self.__number_of_bins_for_chr(hdf5_file, i)
            print("bin for chr: ", bin_for_chr)
            offset_for_chr = self.get_offset_cco_for_chrom(hdf5_file, i)
            j = 0
            while (cur_bin < bin_for_chr):
                start_offset = self.get_offset_for_bin(hdf5_file, cur_bin)
                # borders (cur_bin, cur_bin + batch_size]
                end_offset = self.get_offset_for_bin(hdf5_file, cur_bin + batch_size)
                cur_pixels = pixels[start_offset : min(end_offset, offset_for_chr)]

                cur_pixels = cur_pixels[((cur_pixels["bin1_id"] <= cur_pixels["bin2_id"]) & (
                    cur_pixels["bin1_id"] + (self.up_boader / self.hic_meta.bin_size) >= cur_pixels["bin2_id"]))]
                # print(cur_pixels.values.tolist())

                cur_batch = min(batch_size, bin_for_chr - cur_bin)
                vals, indexies = generate_hic_from_tuples(
                    cur_pixels.values.tolist(), cur_batch, last_size, cur_bin)
                last_size = indexies[len(indexies) - 1] // 2

                self.__save_values(j, i, indexies, vals, chr_name)
                cur_bin += cur_batch
                print(cur_bin)
                j = 1

        print("time: ", time.time() - start_time)

    def get_cmr_sizes(self):
        cmat = cooler.Cooler(f"{self.path}::/resolutions/1000")
        chromsizes = cmat.chromsizes
        print("ajlkfhfdkjhsfj: ", chromsizes[1])

    def get_lines(self, chr: int, start: int, WINDOW_SIZE: int):
        # надо нормально определять пересечения
        chr_name = self.get_chr_name(chr + 1)
        end_pos = math.ceil((start + WINDOW_SIZE) / self.hic_meta.bin_size)
        start_pos = start // self.hic_meta.bin_size
        cur_name = self.vals_name + chr_name + ".bin"
        start_in_file = self.hic_meta.chromsizes[chr].start_pos + start_pos
        # print(start_in_file)
        return read_convert_to_hic_matrix(
            self.indexes_name, cur_name, start_in_file, end_pos - start_pos)

    def get_offset_for_bin(self, group, bin_id : int):
        return group['resolutions'][self.hic_meta.bin_size]['indexes']['bin1_offset'][bin_id]

    def get_offset_cco_for_chrom(self, group, chr : int):
        first_bin1 = self.__number_of_bins_for_chr(group, chr)
        return group['resolutions'][self.hic_meta.bin_size]['indexes']['bin1_offset'][first_bin1] - 1

    def __number_of_bins_for_chr(self, group, chr : int):
        return group['resolutions'][self.hic_meta.bin_size]['indexes']['chrom_offset'][chr]
    
    def debug(self):
        cmat = cooler.Cooler(f"{self.path}::/resolutions/1000")
        chr_name = "chr1"
        chromsizes = cmat.chromsizes

        with h5py.File(self.path, 'r') as file:
            print(file['resolutions']['1000']['indexes']['chrom_offset'][1])
            print("тут мы узнали, какой бин первый во второй хромосоме -> 248957")
            print(file['resolutions']['1000']['indexes']['bin1_offset'][248957])
            print("тут мы узнали, какой отступ у первого бина во второй хромосоме -> 204603218")
            print(file['resolutions']['1000']['pixels']['count'][204603218])
            print("pixels: ", cmat.pixels()[204603218])
            print("тут мы нашли первую запись CCO для второй хромосомы")

    def get_name(self):
        return "hi_c_limit"