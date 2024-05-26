from .Characteristics import Characteristic
import h5py
import threading
import time
import numpy as np
import pandas as pd
from enum import Enum
import os
import math

class HiCCollerFormat():
    """
    Scheme:
    ├── chroms
    │   ├── length (24,) int32
    │   └── name (24,) |S64
    ├── bins
    │   ├── chrom (3088281,) int32
    │   ├── start (3088281,) int32
    │   ├── end (3088281,) int32
    │   └── weight (3088281,) float64
    ├── pixels
    │   ├── bin1_id (271958554,) int64
    │   ├── bin2_id (271958554,) int64
    │   └── count (271958554,) int32
    └── indexes
        ├── bin1_offset (3088282,) int64
        └── chrom_offset (25,) int64
    """
    RESOLUTIONS = 'resolutions/'
    CHROMS_LENGTH = '/chroms/length'
    PIXELS_BIN1_ID = '/pixels/bin1_id'
    PIXELS_BIN2_ID = '/pixels/bin2_id'
    PIXELS_COUNT = '/pixels/count'
    INDEXES_BINS = '/indexes/bin1_offset'
    INDEXES_CH_OFFSET = '/indexes/chrom_offset'
    BINS_CHR = '/bins/chrom'
    BINS_START = '/bins/start'
    BINS_END = '/bins/end'
    BINS_WEIGHT = '/bins/weight'

    def __init__(self, hic_path: str) -> None:
        self.hic_path = hic_path
        self.group = h5py.File(hic_path, 'r')

    def __get_prefix(self) -> str:
        """
        The mcool file has a split into different resolutions.
        """
        if self.hic_path.endswith(".cool"):
            return ''
        elif self.hic_path.endswith(".mcool"):
            return f'{self.RESOLUTIONS}{1_000}'

    def __get_data(self, directory: str, start: np.int64,
                   end: np.int64) -> pd.DataFrame:
        return self.group[f'{self.__get_prefix()}{directory}'][start: end]

    def get_offset(self, chr: int) -> np.int64:
        """
        Offsets for each chromosome corresponds to the first bin_id_1 in the chromosome.
        Len of offset : number of chr + 1

        Parameters
        ----------
        chr: number of chromosome. [0 ... number of chr + 1]
        """
        return self.group[f'{self.__get_prefix()}{self.INDEXES_CH_OFFSET}'][chr]

    def get_chrom_offsets(self) -> pd.core.series.Series:
        """
        All offsets for each chromosome

        """
        return self.group[f'{self.__get_prefix()}{self.INDEXES_CH_OFFSET}']

    def get_offset_for_bin1(self, bin_id: np.int64) -> np.int64:
        """
        Offsets for each bin_id_1. CCO records are sorted by bin_id1.

        """
        return self.group[f'{self.__get_prefix()}{self.INDEXES_BINS}'][bin_id]

    def get_pixels(self, start: np.int64, end: int) -> pd.DataFrame:
        """
        In Coller format pixels represents as:
         pixels
            ├── bin1_id (271958554,) int64
            ├── bin2_id (271958554,) int64
            └── count (271958554,) int32

        Parameters
        ----------
        start: start position in pixels
        """
        bin1_ids = self.__get_data(self.PIXELS_BIN1_ID, start, end)
        bin2_ids = self.__get_data(self.PIXELS_BIN2_ID, start, end)
        counts = self.__get_data(self.PIXELS_COUNT, start, end)
        return pd.DataFrame(
            {'bin1_id': bin1_ids, 'bin2_id': bin2_ids, 'count': counts})
    
    def get_bin1_number(self) -> np.int64:
        return len(self.group[f'{self.__get_prefix()}{self.INDEXES_BINS}'])


    def get_bin_info(self, bin_id_start: np.int64,
                     bin_id_end: np.int64) -> pd.DataFrame:
        """
            bins
            ├── chrom (3088281,) int32
            ├── start (3088281,) int32
            ├── end (3088281,) int32
            └── weight (3088281,) float64
        """
        chroms = self.__get_data(self.BINS_CHR, bin_id_start, bin_id_end)
        starts = self.__get_data(self.BINS_START, bin_id_start, bin_id_end)
        ends = self.__get_data(self.BINS_END, bin_id_start, bin_id_end)
        weights = self.__get_data(self.BINS_WEIGHT, bin_id_start, bin_id_end)
        return pd.DataFrame(
            {'chrom': chroms, 'start': starts, 'end': ends, 'weight': weights})

    def get_bin_size(self) -> int:
        """
        return bin size
        """
        return self.group.attrs['bin-size']
    
class HICFull():
    """
    Scheme:
    ├── meta
    │   ├── square_row, int32
    │   ├── square_column, int32
    │   ├── el_counts, int32
    │   └── start, int64
    └── pixels
    │   ├── row, int64
    │   ├── column,   int64
    │   └── data,  int32
    └── indexes
        └── chrom_offset (25,) int64
    """

    def __init__(self, hic_path: str, bin_size: int) -> None:
        self.hic_path = hic_path
        if os.path.exists(hic_path):
            print(f"File {hic_path} exist.")
            self.group = h5py.File(hic_path, 'r')
        else:
            print(f"File {hic_path} don't find.")
            self.generate_new_file()
        self.bin_size = bin_size
        self.square_size = 300

    def __create_dataset(self, group: h5py._hl.group.Group,
                         d_name: str, d_type: type) -> None:
        group.create_dataset(
            d_name, data=[], chunks=True, maxshape=(
                None,), dtype=d_type)
        print("create dataset")

    def generate_new_file(self) -> None:
        file = h5py.File(self.hic_path, "w", libver='latest', swmr=True)
        pixels_grp = file.create_group("pixels")
        self.__create_dataset(pixels_grp, "bin1_id", np.int64)
        self.__create_dataset(pixels_grp, "bin2_id", np.int64)
        self.__create_dataset(pixels_grp, "data", np.int32)

        meta_grp = file.create_group("meta")
        self.__create_dataset(meta_grp, "square_row", np.int32)
        self.__create_dataset(meta_grp, "square_column", np.int32)
        self.__create_dataset(meta_grp, "el_counts", np.int32)
        self.__create_dataset(meta_grp, "start", np.int64)
        file.swmr_mode = True

        file.close()

    def __update_dataset(self, file: h5py.File, path: str,
                         add_dataset: pd.core.series.Series) -> None:
        dataset = file[path]
        add_len = len(add_dataset)
        dataset.resize((file[path].shape[0] + add_len), axis=0)
        dataset[-add_len:] = add_dataset

        # print("Data stored pixels/start:", dataset[:])

    def update_with_new_square(
            self, square_row: int, square_column: int, start: int, square_size : int,  pixels: pd.DataFrame) -> None:
        print(type(pixels['bin1_id']))
        with h5py.File(self.hic_path, "a") as file:
            ln_elm = len(pixels['count'])
            if (ln_elm > 0):
                self.__update_dataset(file, 'pixels/bin1_id', pixels['bin1_id'].sub(square_row * square_size))
                self.__update_dataset(file, 'pixels/bin2_id', pixels['bin2_id'].sub(square_column * square_size ))
                self.__update_dataset(file, 'pixels/data', pixels['count'])
            self.__update_dataset(
                file,
                'meta/square_row',
                pd.Series(
                    [square_row]))
            self.__update_dataset(
                file,
                'meta/square_column',
                pd.Series(
                    [square_column]))
            self.__update_dataset(file, 'meta/el_counts', pd.Series([ln_elm]))
            self.__update_dataset(file, 'meta/start', pd.Series([start]))

    def get_meta(self, start_r : int, start_c : int, window_size : int) -> pd.DataFrame:
        s_r_start = start_r // self.square_size
        s_r_end = (start_r + window_size) // self.square_size
        s_c_start = start_c // self.square_size
        s_c_end = (start_c + window_size) // self.square_size

        print("s_r_start: ", s_r_start, "s_r_end: ", s_r_end, "s_c_start: ", s_c_start, "s_c_end: ", s_c_end)
        square_rows = self.group['meta/square_row']
        square_columns = self.group['meta/square_column']
        el_counts = self.group['meta/el_counts']
        starts = self.group['meta/start']
        cur_datafr =  pd.DataFrame(
            {'square_rows': square_rows, 'square_columns': square_columns,
              'el_counts': el_counts, 'starts': starts})
        if (s_r_start < s_c_start):
            return cur_datafr[(s_r_start <= cur_datafr['square_rows']) & (s_r_end >= cur_datafr['square_rows']) &
                            (s_c_start <= cur_datafr['square_columns'])  & (s_c_end >= cur_datafr['square_columns'])]
        else:
            return cur_datafr[(s_c_start <= cur_datafr['square_rows']) & (s_c_end >= cur_datafr['square_rows']) &
                            (s_r_start <= cur_datafr['square_columns'])  & (s_r_end >= cur_datafr['square_columns'])]

    
    def get_pixels(self, start : int, end : int):
        row_in_square = self.group['pixels/bin1_id'][start: end]
        column_in_square = self.group['pixels/bin2_id'][start: end]
        data = self.group['pixels/data'][start: end]
        return pd.DataFrame(
            {'bin1_id': row_in_square, 'bin2_id': column_in_square, 'data': data})




class CharacteristicFullHiC(Characteristic):
    def __init__(self, path : str, hic_path : str = ""):
        super().__init__(
            path,
        )
        # first size of squares will be 256
        self.square_size = 300
        self.hic_path = hic_path
        if (self.hic_path == ""):
            self.my_h5 = HICFull(self.path + ".h5", 1_000)

    def preprocess(self):
        # preprocessing in progress
        hic_prev = HiCCollerFormat(hic_path=self.hic_path)
        my_h5 = HICFull(self.path + ".h5", 1_000)
        cur_chr_offset = hic_prev.get_offset(1)
        cur_start = 0
        cur_row = 0
        cur_column = 0
        while (cur_row * self.square_size < cur_chr_offset):
            # Keep only the top part
            pixels = self.__generate_coo_for_square(hic_prev, cur_row, cur_column)
            # print(len(pixels))
            # print(pixels)
            my_h5.update_with_new_square(square_row=cur_row, 
                                         square_column=cur_column, 
                                         start=cur_start, 
                                         square_size=self.square_size, 
                                         pixels=pixels
                                         )
            cur_start += len(pixels)
            cur_column += 1
            if (cur_column * self.square_size >= cur_chr_offset):
                cur_row += 1
                cur_column = cur_row

    def __generate_coo_for_square(
            self, hic_prev: HiCCollerFormat, square_row: int, square_column: int) -> pd.DataFrame:
        """
        Finds all the elements for the square under the number.
        +1 because fun shows the indentation for the current bin_id

        Parameters
        ----------
        square: number of square
        """
        print("square row: ", square_row, " square column: ", square_column)
        start = hic_prev.get_offset_for_bin1(square_row * self.square_size)
        end = hic_prev.get_offset_for_bin1((square_row + 1) * self.square_size + 1)
        pixels = hic_prev.get_pixels(start, end)
        return pixels[(square_column * self.square_size <= pixels["bin2_id"])
                      & (pixels["bin2_id"] <= (square_column + 1) * self.square_size)]

    def get_lines(self, chr: int, start_0 : int, start_1 : int, window_size : int):
        new_start_0 = start_0 // self.my_h5.bin_size
        new_start_1 = start_1 // self.my_h5.bin_size
        new_window_size = math.ceil((start_0 + window_size) / self.my_h5.bin_size) - new_start_0

        print("new start 0:", new_start_0, "new start 1:", new_start_1)
        squares = self.my_h5.get_meta(new_start_0, new_start_1, new_window_size)
        matrix = np.zeros((new_window_size, new_window_size), dtype=int)
        print(squares)
        for _, square in squares.iterrows():
            start = square['starts']
            count_of_els = square['el_counts']
            # print(start, count_of_els)
            cur_pixels = self.my_h5.get_pixels(start, start + count_of_els)

            cur_pixels['c_r_1'] = cur_pixels['bin1_id'] +  square["square_rows"] * self.my_h5.square_size - new_start_0
            cur_pixels['c_c_1'] = cur_pixels['bin2_id'] + square["square_columns"] * self.my_h5.square_size - new_start_1 

            cur_pixels['c_r_2'] = cur_pixels['bin2_id'] + square["square_columns"] * self.my_h5.square_size - new_start_0
            cur_pixels['c_c_2'] = cur_pixels['bin1_id'] +  square["square_rows"] * self.my_h5.square_size - new_start_1 
            
            cur_pixels_1 = cur_pixels[((cur_pixels['c_r_1'] < new_window_size) &
                                    (cur_pixels['c_c_1']< new_window_size) &
                                    (0 <= cur_pixels['c_r_1']) &
                                    (0 <= cur_pixels['c_c_1']))]
            
            
            cur_pixels_2 = cur_pixels[(cur_pixels['c_r_2'] < new_window_size) &
                                     (cur_pixels['c_c_2'] < new_window_size) &
                                     (0 <= cur_pixels['c_r_2']) &
                                     (0 <= cur_pixels['c_c_2'])] 
            
            # del cur_pixels

            matrix[cur_pixels_1['c_r_1'], cur_pixels_1['c_c_1']] = cur_pixels_1['data']
            matrix[cur_pixels_2['c_r_2'], cur_pixels_2['c_c_2']] = cur_pixels_2['data']
            

        print(new_window_size)
        return matrix

