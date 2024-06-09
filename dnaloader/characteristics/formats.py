import h5py
import numpy as np
import pandas as pd

class CoolerFormat():
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
    BINS = 'bins'
    BINS_CHR = '/' + BINS + '/chrom'
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
            return f'{self.RESOLUTIONS}{self.get_bin_size()}'

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
    
    def get_chrom_lengths(self) -> pd.core.series.Series:
        """
        Length of all chromosomes, measured in nucleotides

        """
        return self.group[f'{self.__get_prefix()}{self.CHROMS_LENGTH}']

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
        If the file is of the .mcool type with different resolutions, 
        then function returns the minimum. If.cool - bin_size.
        """
        if self.hic_path.endswith(".cool"):
            self.group.attrs['bin-size']
        elif self.hic_path.endswith(".mcool"):
            return min(self.group['resolutions'])