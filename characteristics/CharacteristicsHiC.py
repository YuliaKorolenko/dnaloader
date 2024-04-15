from .Characteristics import Characteristic
import cooler
import time
import numpy as np
import math
from characteristics.bwhelper import generate_array_from_tuples, write_array_to_file, read_numbers_from_file, read_convert_to_hic_matrix, generate_hic_from_tuples
import os
from dataclasses import dataclass
import pickle


@dataclass
class HiC_Meta:
    bin_size: int


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
        self.vals_name = os.path.join(self.folder_res, "vals.bin")
        self.meta_name = os.path.join(self.folder_res, "hic_meta.pickle")
        with open(self.meta_name, 'rb') as file:
            self.hic_meta = pickle.load(file)
        self.up_boader = 1_000_000

    def preprocess(self):
        if not os.path.exists(self.folder_res):
            os.mkdir(self.folder_res)
        cmat = cooler.Cooler(f"{self.path}::/resolutions/1000")

        self.hic_meta = HiC_Meta(bin_size=cmat.info['bin-size'])
        with open(self.meta_name, 'wb') as file:
            pickle.dump(self.hic_meta, file)

        chromsizes = cmat.chromsizes

        pixels = cmat.pixels()

        start_time = time.time()
        batch_size = 100_000_000
        j = 0
        last_size = 0
        last_max_row = 0

        prev_row = 0
        while (j < chromsizes["chr1"]):
            cur_pixels = pixels[j: min(j + batch_size, chromsizes["chr1"])]
            max_row = cur_pixels['bin1_id'].max()
            if (j + len(cur_pixels) < chromsizes["chr1"]):
                cur_pixels = cur_pixels[(cur_pixels["bin1_id"] < max_row)]
            len_of_pixels = len(cur_pixels)
            cur_pixels = cur_pixels[((cur_pixels["bin1_id"] <= cur_pixels["bin2_id"]) & (
                cur_pixels["bin1_id"] + (self.up_boader / self.hic_meta.bin_size) >= cur_pixels["bin2_id"]))]

            vals, indexies = generate_hic_from_tuples(
                cur_pixels.values.tolist(), max_row - last_max_row, last_size, last_max_row)
            last_size = indexies[len(indexies) - 1] // 2
            last_max_row = max_row

            mode = b'a'
            if (j == 0):
                mode = b'wb'

            if (j == 0):
                write_array_to_file(
                    self.indexes_name, np.array(
                        indexies, dtype=np.int32), mode)
            else:
                write_array_to_file(self.indexes_name, np.array(
                    indexies[1:], dtype=np.int32), mode)

            write_array_to_file(
                self.vals_name, np.array(
                    vals, dtype=np.int32), mode)
            del indexies
            del vals

            print(j)
            j += len_of_pixels

        print("time: ", time.time() - start_time)

    def get_lines(self, chr: int, start: int, WINDOW_SIZE: int):
        # надо нормально определять пересечения
        end_pos = math.ceil((start + WINDOW_SIZE) / self.hic_meta.bin_size)
        print("end pos: ", end_pos)
        start_pos = start // self.hic_meta.bin_size
        print("start pos: ", start_pos)
        return read_convert_to_hic_matrix(
            self.indexes_name, self.vals_name, start_pos, end_pos - start_pos)

    def get_name(self):
        return "hi_c_limit"
