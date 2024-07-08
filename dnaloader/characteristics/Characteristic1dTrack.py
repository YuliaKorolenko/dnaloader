import os
import h5py
import numpy as np
import pyBigWig
import math
import numpy as np
from tqdm import tqdm
import psutil
from dnaloader.characteristics.bwhelper import generate_array_from_tuples

class CharacteristicPreprocess():

    def __init__(self, prev_path: str,  res_path: str,):
        self.prev_path = prev_path
        self.res_path = res_path

    def get_chr_name(self, i: int):
        chr = "chr"
        if i == 23:
            chr += "X"
        elif i == 24:
            chr += "Y"
        else:
            chr += str(i)
        return (chr)

class Format1D():
    """
    Scheme:
    ├── pixels
    │   ├── offsets, int64
    │   └── data,  int32
    └── indexes
        └── chroms_offset (25,) int32
    """
    PIX = "pixels"
    PIX_OFFSETS = PIX + "/offsets"
    PIX_DATA = PIX + "/data"

    IND = "indexes"
    IND_OFFSET = IND + "/chroms_offset"

    def __init__(self, hic_path: str) -> None:
        self.hic_path = hic_path

    def __create_dataset(self, group: h5py._hl.group.Group,
                         d_name: str, d_type: type) -> None:
        """
        Create dataset in the group.

        Parameters
        ----------
        group: the group where the dataset is created
        d_name: name of dataset
        d_type: types of values in the dataset
        """
        group.create_dataset(
            d_name, data=[], chunks=True, maxshape=(
                None,), dtype=d_type)
        
    def __get_short(self, fullname : str) -> str:
        """
        Returns the name of the dataset, without the name of the group.

        Parameters
        ----------
        fullname: name of dataset with group name
        """
        return fullname.split("/")[-1]

    def generate_new_file(self) -> None:
        """
        Generate file with scheme.
        """
        self.group = h5py.File(self.hic_path, "w", libver='latest', swmr=True)
        pixels_grp = self.group.create_group(self.PIX)
        self.__create_dataset(pixels_grp, self.__get_short(self.PIX_OFFSETS), np.int64)
        self.__create_dataset(pixels_grp, self.__get_short(self.PIX_DATA), np.int32)

        ind_grp = self.group.create_group(self.IND)
        self.__create_dataset(ind_grp, self.__get_short(self.IND_OFFSET), np.int32)
        self.group.swmr_mode = True

        self.group.close()

    def __update_dataset(self, dataset: h5py.Dataset,
                         add_dataset) -> None:
        """
        Adding values to an dataset
        
        Parameters
        ----------
        dataset: the dataset that needs to be increased
        add_dataset: the array that is appended to the dataset
        """
        add_len = len(add_dataset)
        dataset.resize((dataset.shape[0] + add_len), axis=0)
        dataset[-add_len:] = add_dataset

    def update_with_new_part(self, pixel_offsets, pixel_data) -> None:
        with h5py.File(self.hic_path, "a") as file:
            self.__update_dataset(file[self.PIX_OFFSETS], pixel_offsets)
            if (len(pixel_data) > 0):
                self.__update_dataset(file[self.PIX_DATA], pixel_data)

        del pixel_offsets, pixel_data

    def update_meta(self, ind_offsets, bin_size: int):
        with h5py.File(self.hic_path, "a") as file:
            file.attrs['bin_size'] = 1
            self.__update_dataset(file[self.IND_OFFSET], ind_offsets)
        

    
class Characteristic1dPreprocess(CharacteristicPreprocess):
    def __init__(self, prev_path: str, res_path: str = ""):
        super().__init__(
            prev_path=prev_path,
            res_path=res_path
        )

    def preprocess(self):
        self.create_new_file()

        # generate all bw paths
        file_paths = []
        with open(self.prev_path, 'r') as file:
            for line in file:
                file_paths.append(line.strip())

        bw = pyBigWig.open(file_paths[0])
        new_format = self.__generate_offset_new_format(bw.chroms())
        self.format1d.update_meta(new_format, 1)
        self.preprocess_subsequence(file_paths, new_format)

    def __generate_offset_new_format(self, bw_chroms):
        new_format_offset = []
        for i in range(1, 24 + 1):
            chr = "chr"
            if i == 23:
                chr += "X"
            elif i == 24:
                chr += "Y"
            else:
                chr += str(i)

            chrom_size = bw_chroms[chr]
            new_format_offset.append(chrom_size)
        return new_format_offset

    def save_information(self, arr_col_data, chr_name: str,
                         cur_pos: int, index_in_row_cur, chr_num: int):
        self.format1d.update_with_new_part(index_in_row_cur, arr_col_data)
        del index_in_row_cur
        del arr_col_data
    
    def get_max_array_size(self):
        available_memory = psutil.virtual_memory().available
        dtype_size = np.dtype(np.int64).itemsize
        max_elements = available_memory // dtype_size
        # add bw count
        return max_elements // 6 // 10

    def preprocess_subsequence(self, file_paths, ind_offsets):
        for num_chr in range(1, 24 + 1):
            last_post_in_file = 0
            cur_pos = 0
            chrom_size = ind_offsets[num_chr - 1]
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

    def create_new_file(self):
        self.format1d = Format1D(self.res_path)
        self.format1d.generate_new_file()
        