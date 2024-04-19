from .Characteristics import Characteristic
import h5py
import threading
import time
import numpy as np
import pandas as pd
from characteristics.bwhelper import read_numbers_from_file
from enum import Enum
import os


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

    def __init__(self, hic_path: str, bin_size: int) -> None:
        self.hic_path = hic_path
        self.group = h5py.File(hic_path, 'r')
        self.bin_size = bin_size

    def __get_prefix(self) -> str:
        """
        The mcool file has a split into different resolutions.
        """
        if self.hic_path.endswith(".cool"):
            return ''
        elif self.hic_path.endswith(".mcool"):
            return f'{self.RESOLUTIONS}{self.bin_size}'

    def __get_data(self, directory: str, start: np.int64,
                   end: np.int64) -> pd.DataFrame:
        return self.group[f'{self.__get_prefix()}{directory}'][start: end]

    def get_offset(self, chr: int) -> np.int64:
        """
        Offsets for each chromosome corresponds to the first bin_id_1 in the chromosome.
        Len of offset : number of chr + 1 
        :param chr: number of chromosome. [0 ... number of chr + 1]
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

        :param start: start position in pixels
        """
        bin1_ids = self.__get_data(self.PIXELS_BIN1_ID, start, end)
        bin2_ids = self.__get_data(self.PIXELS_BIN2_ID, start, end)
        counts = self.__get_data(self.PIXELS_COUNT, start, end)
        return pd.DataFrame(
            {'bin1_id': bin1_ids, 'bin2_id': bin2_ids, 'count': counts})

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
    │   ├── start, int64
    │   ├── end,   int64
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

    def __create_dataset(self, group : h5py._hl.group.Group, d_name : str, d_type : type) -> None:
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

        print("Data stored pixels/start:", dataset[:])

    def update_with_new_square(
            self, square_row: int, square_column: int,  start: int, pixels: pd.DataFrame) -> None:
        print(type(pixels['bin1_id']))
        with h5py.File(self.hic_path, "a") as file:
            ln_elm = len(pixels['count'])
            self.__update_dataset(file, 'pixels/bin1_id', pixels['bin1_id'])
            self.__update_dataset(file, 'pixels/bin2_id', pixels['bin2_id'])
            self.__update_dataset(file, 'pixels/data', pixels['count'])
            self.__update_dataset(file, 'meta/square_row', pd.Series([square_row]))
            self.__update_dataset(file, 'meta/square_column', pd.Series([square_column]))
            self.__update_dataset(file, 'meta/el_counts', pd.Series([ln_elm]))
            self.__update_dataset(file, 'meta/start', pd.Series([start]))


class CharacteristicFullHiC(Characteristic):
    def __init__(self, hic_path: str):
        super().__init__(
            hic_path,
        )

    def preprocess(self):
        # preprocessing in progress
        hic_prev = HiCCollerFormat(hic_path=self.path, bin_size=1_000)
        pixels = self.__generate_coo_for_square(hic_prev, 0, 0)
        my_h5 = HICFull("loader.h5", 1_000)
        my_h5.update_with_new_square(0, 0, 0, pixels)

    def __generate_coo_for_square(
            self, hic_prev: HiCCollerFormat, square_row: int, square_column : int) -> pd.DataFrame:
        """
        Finds all the elements for the square under the number.
        +1 because fun shows the indentation for the current bin_id

        :param square: number of square
        """
        # first size of squares will be 256
        size = 256
        start = hic_prev.get_offset_for_bin1(square_row * size)
        end = hic_prev.get_offset_for_bin1((square_row + 1) * size + 1)
        pixels = hic_prev.get_pixels(start, end)
        return pixels[(square_column * size <= pixels["bin2_id"]) & (pixels["bin2_id"] <= (square_column + 1) * size)]

    def get_lines(self):
        print(self.group['pixels/bin1_id'][0:100])
