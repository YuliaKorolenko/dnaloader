from .Characteristics import Characteristic
from .formats import CoolerFormat
import h5py
import numpy as np
import pandas as pd
import os
import math

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
        ├── chroms_offset (25,) int32
        └── chroms_length (25,) int32
    """
    META = "meta"
    META_S_ROW = META + "/square_row"
    META_S_COL = META + "/square_column"
    META_EL_COUNT = META +"/el_counts"
    META_START = META + "/start"

    PIX = "pixels"
    PIX_BIN1 = PIX + "/bin1_id"
    PIX_BIN2 = PIX + "/bin2_id"
    PIX_DATA = PIX + "/data"

    IND = "indexes"
    IND_OFFSET = IND + "/chroms_offset"
    IND_LENGTH = IND + "/chroms_length"

    def __get_short(self, fullname : str):
        return fullname.split("/")[-1]

    def __init__(self, hic_path: str) -> None:
        self.hic_path = hic_path
        if os.path.exists(hic_path):
            print(f"File {hic_path} exist.")
            self.group = h5py.File(hic_path, 'r')
        else:
            print(f"File {hic_path} don't find.")
            self.__generate_new_file()


    def __create_dataset(self, group: h5py._hl.group.Group,
                         d_name: str, d_type: type) -> None:
        group.create_dataset(
            d_name, data=[], chunks=True, maxshape=(
                None,), dtype=d_type)

    def __generate_new_file(self) -> None:
        self.group = h5py.File(self.hic_path, "w", libver='latest', swmr=True)
        pixels_grp = self.group.create_group(self.PIX)
        print("aaaa ", self.__get_short(self.PIX_BIN1))
        self.__create_dataset(pixels_grp, self.__get_short(self.PIX_BIN1), np.int64)
        self.__create_dataset(pixels_grp, self.__get_short(self.PIX_BIN2), np.int64)
        self.__create_dataset(pixels_grp, self.__get_short(self.PIX_DATA), np.int32)

        meta_grp = self.group.create_group(self.META)
        self.__create_dataset(meta_grp, self.__get_short(self.META_S_ROW), np.int32)
        self.__create_dataset(meta_grp, self.__get_short(self.META_S_COL), np.int32)
        self.__create_dataset(meta_grp, self.__get_short(self.META_EL_COUNT), np.int32)
        self.__create_dataset(meta_grp, self.__get_short(self.META_START), np.int64)

        print("before creating ind")
        ind_grp = self.group.create_group(self.IND)
        self.__create_dataset(ind_grp, self.__get_short(self.IND_OFFSET), np.int32)
        self.__create_dataset(ind_grp, self.__get_short(self.IND_LENGTH), np.int32)
        self.group.swmr_mode = True

        self.group.close()

    def __update_dataset(self, file: h5py.File, path: str,
                         add_dataset: pd.core.series.Series) -> None:
        dataset = file[path]
        add_len = len(add_dataset)
        dataset.resize((file[path].shape[0] + add_len), axis=0)
        dataset[-add_len:] = add_dataset

    def generate_offsets_and_header(self, hic_prev: CoolerFormat) -> None:
        with h5py.File(self.hic_path, "a") as file:
            self.__update_dataset(file, self.IND_OFFSET, hic_prev.get_chrom_offsets())
            self.__update_dataset(file, self.IND_LENGTH, hic_prev.get_chrom_lengths())
            file.attrs['bin_size'] = hic_prev.get_bin_size()
            file.file.attrs['square_size'] = 300

    def update_with_new_square(
            self, square_row: int, square_column: int, start: int, pixels: pd.DataFrame) -> None:
        print(type(pixels['bin1_id']))
        with h5py.File(self.hic_path, "a") as file:
            ln_elm = len(pixels['count'])
            if (ln_elm > 0):
                self.__update_dataset(file, self.PIX_BIN1, pixels['bin1_id'].sub(square_row * self.get_square_size()))
                self.__update_dataset(file, self.PIX_BIN2, pixels['bin2_id'].sub(square_column * self.get_square_size() ))
                self.__update_dataset(file, self.PIX_DATA, pixels['count'])
            self.__update_dataset(file, self.META_S_ROW, pd.Series([square_row]))
            self.__update_dataset(file, self.META_S_COL,pd.Series([square_column]))
            self.__update_dataset(file, self.META_EL_COUNT, pd.Series([ln_elm]))
            self.__update_dataset(file, self.META_START, pd.Series([start]))

    def get_meta(self, start_r: int, start_c: int,
                 window_size: int) -> pd.DataFrame:
        
        s_r_start = start_r // self.get_square_size()
        s_r_end = math.ceil((start_r + window_size) / self.get_square_size())
        s_c_start = start_c // self.get_square_size()
        s_c_end = math.ceil((start_c + window_size) / self.get_square_size())

        square_rows = self.group[self.META_S_ROW]
        square_columns = self.group[self.META_S_COL]
        el_counts = self.group[self.META_EL_COUNT]
        starts = self.group[self.META_START]
        cur_datafr = pd.DataFrame(
            {'square_rows': square_rows, 'square_columns': square_columns,
             'el_counts': el_counts, 'starts': starts})
        if (s_r_start < s_c_start):
            return cur_datafr[(s_r_start <= cur_datafr['square_rows']) & (s_r_end >= cur_datafr['square_rows']) &
                              (s_c_start <= cur_datafr['square_columns']) & (s_c_end >= cur_datafr['square_columns'])]
        else:
            return cur_datafr[(s_c_start <= cur_datafr['square_rows']) & (s_c_end >= cur_datafr['square_rows']) &
                              (s_r_start <= cur_datafr['square_columns']) & (s_r_end >= cur_datafr['square_columns'])]

    def get_pixels(self, start: int, end: int):
        """
        Returns arays from group pixels.

        Parameters
        ----------
        start: start position in arrays from pixel group
        end: end position in arrays from pixel group
        """
        row_in_square = self.group[self.PIX_BIN1][start: end]
        column_in_square = self.group[self.PIX_BIN2][start: end]
        data = self.group[self.PIX_DATA][start: end]
        return pd.DataFrame(
            {'bin1_id': row_in_square, 'bin2_id': column_in_square, 'data': data})

    def get_chr_length(self, chr: int) -> np.int64:
        """
        Returns length for chromosome. 
        This length is equal to the number of nucleotids in this chromosome.
        Equals to CHROMS_LENGTH in cooler format.

        Parameters
        ----------
        chr: chromosome number [0 ... 24]
        """
        return self.group[self.IND_LENGTH][chr]
    
    def get_chr_offset(self, chr : int) -> np.int64:
        """
        Returns offsets for chromosome. 
        This offset for one chromosome is equal to the total number 
        of non-zero elements in all previous chromosomes
        Equals to INDEXES_CH_OFFSET in cooler format.

        Parameters
        ----------
        chr: chromosome number [0 ... 24]
        """
        return self.group[self.IND_OFFSET][chr]
    
    def get_bin_size(self) -> np.int64:
        return int(self.group.attrs['bin_size'])
    
    def get_square_size(self) -> np.int64:
        # return self.group.attrs['square_size']
        return 300


class CharacteristicFullHiC(Characteristic):
    def __init__(self, path: str, hic_path: str = ""):
        super().__init__(
            path,
        )
        self.hic_path = hic_path
        if (self.hic_path == ""):
            self.full_hic = HICFull(self.path)

    def preprocess(self):
        # preprocessing in progress
        hic_prev = CoolerFormat(hic_path=self.hic_path)

        full_hic = HICFull(self.path)

        hic_prev.get_bin_size()

        full_hic.generate_offsets_and_header(hic_prev=hic_prev)
        
        cur_chr_offset = hic_prev.get_offset(1)
        cur_start = 0
        cur_row = 0
        cur_column = 0
        while (cur_row * full_hic.get_square_size() < cur_chr_offset):
            # Keep only the top part
            pixels = self.__generate_coo_for_square(
                hic_prev, cur_row, full_hic.get_square_size(), cur_column)
            full_hic.update_with_new_square(square_row=cur_row,
                                         square_column=cur_column,
                                         start=cur_start,
                                         pixels=pixels
                                         )
            cur_start += len(pixels)
            cur_column += 1
            if (cur_column * full_hic.get_square_size() >= cur_chr_offset):
                cur_row += 1
                cur_column = cur_row
        self.group.close()
        self.group = h5py.File(self.hic_path, 'r')

    def get_chr_lenght(self, chr: int) -> np.int64:
        return self.full_hic.get_chr_length(chr)

    def __generate_coo_for_square(
            self, hic_prev: CoolerFormat, square_row: int, square_size : int, square_column: int) -> pd.DataFrame:
        """
        Finds all the elements for the square under the number.
        +1 because fun shows the indentation for the current bin_id

        Parameters
        ----------
        square: number of square
        """
        print("square row: ", square_row, " square column: ", square_column)
        start = hic_prev.get_offset_for_bin1(square_row * square_size)
        end = hic_prev.get_offset_for_bin1((square_row + 1) * square_size + 1)
        pixels = hic_prev.get_pixels(start, end)
        return pixels[(square_column * square_size <= pixels["bin2_id"])
                      & (pixels["bin2_id"] <= (square_column + 1) * square_size)]

    def get_lines(self, chr: int, start_0: int,
                  start_1: int, window_size: int) -> np.array:
        start_chr = self.full_hic.get_chr_offset(chr)
        new_start_0 = start_chr + start_0 // self.full_hic.get_bin_size()
        new_start_1 = start_chr + start_1 // self.full_hic.get_bin_size()
        new_window_size = math.ceil(
            (start_chr + start_0 + window_size) / self.full_hic.get_bin_size()) - new_start_0

        squares = self.full_hic.get_meta(
            new_start_0, new_start_1, new_window_size)
        matrix = np.zeros((new_window_size, new_window_size), dtype=int)
        for _, square in squares.iterrows():
            start = square['starts']
            count_of_els = square['el_counts']
            cur_pixels = self.full_hic.get_pixels(start, start + count_of_els)

            cur_pixels['c_r_1'] = cur_pixels['bin1_id'] + \
                square["square_rows"] * self.full_hic.get_square_size() - new_start_0
            cur_pixels['c_c_1'] = cur_pixels['bin2_id'] + \
                square["square_columns"] * self.full_hic.get_square_size()- new_start_1

            cur_pixels['c_r_2'] = cur_pixels['bin2_id'] + \
                square["square_columns"] * self.full_hic.get_square_size()- new_start_0
            cur_pixels['c_c_2'] = cur_pixels['bin1_id'] + \
                square["square_rows"] * self.full_hic.get_square_size() - new_start_1

            cur_pixels_1 = cur_pixels[((cur_pixels['c_r_1'] < new_window_size) &
                                       (cur_pixels['c_c_1'] < new_window_size) &
                                       (0 <= cur_pixels['c_r_1']) &
                                       (0 <= cur_pixels['c_c_1']))]

            cur_pixels_2 = cur_pixels[(cur_pixels['c_r_2'] < new_window_size) &
                                      (cur_pixels['c_c_2'] < new_window_size) &
                                      (0 <= cur_pixels['c_r_2']) &
                                      (0 <= cur_pixels['c_c_2'])]

            del cur_pixels

            matrix[cur_pixels_1['c_r_1'],
                   cur_pixels_1['c_c_1']] = cur_pixels_1['data']
            matrix[cur_pixels_2['c_r_2'],
                   cur_pixels_2['c_c_2']] = cur_pixels_2['data']

        return matrix
