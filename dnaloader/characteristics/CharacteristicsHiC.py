from .Characteristics import Characteristic
import cooler
import time
import numpy as np
import math
from dnaloader.characteristics.bwhelper import write_array_to_file, read_numbers_from_file, read_convert_to_hic_matrix, generate_hic_from_tuples
import os
from dataclasses import dataclass
import pickle
import pandas as pd
from typing import List
from dnaloader.common import Chromosome_Info
import h5py
from pympler import muppy, summary
from .CharacteristicFullHiC import HiCCollerFormat


@dataclass
class HiC_Meta:
    bin_size: int
    chromsizes: List[Chromosome_Info]


class CharacteristicHiCColer(Characteristic):
    def __init__(self, hic_path: str, bin_size : int):
        super().__init__(
            hic_path,
        )
        if hic_path.endswith(".mcool"):
            self.c_matrix = cooler.Cooler(f"{hic_path}::/resolutions/{bin_size}")
        elif hic_path.endswith(".cool"):
            self.c_matrix = cooler.Cooler(f"{hic_path}")
        else:
            raise ValueError("Wrong path, only files are supported with .cool and .mcool extensions")

    def get_lines(self, chr: int, start: int, WINDOW_SIZE: int):
        return self.c_matrix.matrix(balance=False).fetch(
            (self.get_chr_name(chr + 1), start, start + WINDOW_SIZE))

    def get_line_d(self, chr: int, start_1: int, start_2 : int, WINDOW_SIZE: int):
        return self.c_matrix.matrix(balance=False).fetch(
            (self.get_chr_name(chr + 1), start_1, start_1 + WINDOW_SIZE),
            (self.get_chr_name(chr + 1), start_2, start_2 + WINDOW_SIZE)
            )


    def get_name(self):
        return "hi_c"


class Characteristic2dWithLimit(Characteristic):
    def __init__(self, folder_res: str, loader_type: str, hic_path="", bin_size=1_000):
        super().__init__(
            folder_res,
        )
        self.folder_res = folder_res
        self.indexes_name = os.path.join(self.folder_res, "indexies")
        self.vals_name = os.path.join(self.folder_res, "vals")
        self.meta_name = os.path.join(self.folder_res, "hic_meta.pickle")
        self.chrom_count = 24
        self.index_in_row = np.empty((0), dtype=np.int64)
        self.type = loader_type
        self.hic_path = hic_path
        self.up_boader = 1_000_000
        if (hic_path == ""):
            with open(self.meta_name, 'rb') as file:
                self.hic_meta = pickle.load(file)
            if (self.type == "medium"):
                self.index_in_row = np.load(self.indexes_name + '.npy')
        else:
            self.preprocess(bin_size)

    def __preprocess_meta(
            self, chromsizes: pd.core.series.Series, bin_size: int) -> None:
        """
        This is where meta information is created for hic with limit.
        + i since there is 1 more element in the indexies array for each chromosome than the size

        Parameters
        ----------
        chromsizes: offsets for bins
        bin_size: resolution
        """
        chrom_lists = []
        for i in range(0, self.chrom_count):
            chrom_lists.append(
                Chromosome_Info(
                    number=i,
                    start_pos=chromsizes[i] + i,
                    lenght=chromsizes[i + 1] - chromsizes[i]))
        self.hic_meta = HiC_Meta(bin_size=bin_size, chromsizes=chrom_lists)
        with open(self.meta_name, 'wb') as file:
            pickle.dump(self.hic_meta, file)

    def __save_values(self, is_start: bool, i: int, indexies,
                      vals, chr_name: str) -> None:
        mode = b'a'
        if (is_start):
            mode = b'wb'

        if (self.type == "hard"):
            if (is_start):
                write_array_to_file(
                    self.indexes_name + chr_name + ".bin", np.array(
                        indexies, dtype=np.int32), mode)
            else:
                write_array_to_file(self.indexes_name + chr_name + ".bin", np.array(
                    indexies[1:], dtype=np.int32), mode)
        else:
            self.index_in_row = np.append(
                self.index_in_row[:-1], indexies)
            np.save(self.indexes_name + '.npy', self.index_in_row)

        write_array_to_file(
            self.vals_name + chr_name + ".bin", np.array(
                vals, dtype=np.int32), mode)
        del indexies
        del vals, is_start, chr_name, mode

    def check_memory(self):
        all_objects = muppy.get_objects()
        sum1 = summary.summarize(all_objects)
        # Prints out a summary of the large objects
        summary.print_(sum1)
        # Get references to certain types of objects such as dataframe
        dataframes = [ao for ao in all_objects if isinstance(ao, pd.DataFrame)]
        for d in dataframes:
            print(d.columns.values)
            print(len(d))

    def preprocess(self, bin_size : int) -> None:
        """
        A two-dimensional track around the main diagonal is preprocessed

        """
        if not os.path.exists(self.folder_res):
            os.mkdir(self.folder_res)

        coller_format = HiCCollerFormat(hic_path=self.hic_path)
        # coller_format.get_bin_size()
        self.__preprocess_meta(coller_format.get_chrom_offsets(), bin_size)

        start_time = time.time()
        # TODO: calculate batch size using free ram size
        batch_size = 20_000
        cur_bin = coller_format.get_offset(1 - 1)

        for i in range(1, self.chrom_count + 1):
            chr_name = self.get_chr_name(i)
            last_size = 0
            print(chr_name)
            # requesting offset for the next chromosome
            is_start_new_chr = True
            bin_for_chr = coller_format.get_offset(i)
            print("bin for chr: ", bin_for_chr)
            offset_for_chr = coller_format.get_offset_for_bin1(bin_for_chr) - 1
            while (cur_bin < bin_for_chr):
                start_offset = coller_format.get_offset_for_bin1(cur_bin)
                # borders (cur_bin, cur_bin + batch_size]
                end_offset = coller_format.get_offset_for_bin1(
                    min(cur_bin + batch_size, coller_format.get_bin1_number()) - 1)
                cur_pixels = coller_format.get_pixels(
                    start_offset, min(end_offset, offset_for_chr))

                cur_pixels = cur_pixels[((cur_pixels["bin1_id"] <= cur_pixels["bin2_id"]) & (
                    cur_pixels["bin1_id"] + (self.up_boader / self.hic_meta.bin_size) >= cur_pixels["bin2_id"]))]

                cur_batch = min(batch_size, bin_for_chr - cur_bin)
                vals, indexies = generate_hic_from_tuples(
                    cur_pixels.values.tolist(), cur_batch, last_size, cur_bin)
                del cur_pixels
                last_size = indexies[len(indexies) - 1] // 2

                self.__save_values(
                    is_start_new_chr, i, indexies, vals, chr_name)
                del indexies
                del vals
                cur_bin += cur_batch
                print(cur_bin)
                is_start_new_chr = False

        print("time: ", time.time() - start_time)

    def get_lines(self, chr: int, start: int, WINDOW_SIZE: int):
        """
        Returns an submatrix aroun main diagonal
        on request: chromosome, starting position and window size

        Parameters
        ----------
        chr: the number of the chromosome from which the submatrix was requested
        start: start position in chromosome in nucleotids
        WINDOW_SIZE: size of submatrix in nucleotids
        """
        chr_name = self.get_chr_name(chr + 1)
        end_pos = math.ceil((start + WINDOW_SIZE) / self.hic_meta.bin_size)
        start_pos = start // self.hic_meta.bin_size
        cur_name = self.vals_name + chr_name + ".bin"
        start_in_file = self.hic_meta.chromsizes[chr].start_pos + start_pos
        end_in_file = self.hic_meta.chromsizes[chr].start_pos + end_pos
        if (self.type == "hard"):
            indixies = read_numbers_from_file(
                self.indexes_name + chr_name + ".bin",
                end_pos - start_pos + 1,
                start_pos)
        else:
            indixies = self.index_in_row[start_in_file: end_in_file]
        return read_convert_to_hic_matrix(
            cur_name, np.array(indixies), end_pos - start_pos)
    
    def get_chr_lenght(self, chr : int):
        return self.hic_meta.chromsizes[chr].lenght

    def get_name(self):
        return "hi_c_limit"
