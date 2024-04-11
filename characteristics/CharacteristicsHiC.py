from .Characteristics import Сharacteristic
import cooler
import time
import numpy as np
from characteristics.bwhelper import generate_array_from_tuples, write_array_to_file, read_numbers_from_file, read_convert_to_hic_matrix, generate_hic_from_tuples
import os

class СharacteristicHiCColer(Сharacteristic):
    def __init__(self, hic_path : str):
        super().__init__(
            hic_path,
        )
        self.c_matrix = cooler.Cooler(f"{hic_path}::/resolutions/1000")

    def get_lines(self, chr : int, start : int, WINDOW_SIZE : int):
        return self.c_matrix.matrix(balance=True).fetch((self.get_chr_name(chr + 1), start, start + WINDOW_SIZE))

    def get_name(self):
        return "hi_c"

class CharacteristicHiCWithLimit(Сharacteristic):
    def __init__(self, hic_path : str):
        super().__init__(
            hic_path,
        )
    
    def preprocess(self):
        if not os.path.exists('lalal'):
            os.mkdir('lalal')
        hic_file = "/home/ojpochemy/SamplerBigWig/hi_c/4DNFI9FVHJZQ.mcool"
        cmat = cooler.Cooler(f"{hic_file}::/resolutions/1000")

        chromsizes = cmat.chromsizes

        pixels = cmat.pixels()

        start_time = time.time() 
        batch_size = 50_000_000
        j = 0
        last_size = 0
        last_max_row = 0

        prev_row = 0
        while (j < chromsizes["chr1"]):
            cur_pixels = pixels[j : min(j + batch_size, chromsizes["chr1"])]
            max_row = cur_pixels['bin1_id'].max()
            if (j + len(cur_pixels) < chromsizes["chr1"]):
                cur_pixels = cur_pixels[(cur_pixels["bin1_id"] < max_row)]
            len_of_pixels = len(cur_pixels)
            cur_pixels = cur_pixels[((cur_pixels["bin1_id"] <= cur_pixels["bin2_id"])&(cur_pixels["bin1_id"] + (1_000_000 / 1_000) >= cur_pixels["bin2_id"]))]
            
            vals, indexies = generate_hic_from_tuples(cur_pixels.values.tolist(), max_row - last_max_row, last_size, last_max_row)
            last_size = indexies[len(indexies) - 1] // 2
            last_max_row = max_row

            mode =  b'a'
            if (j == 0):
                mode = b'wb'

            if (j == 0):
                write_array_to_file("lalal/indexies.bin", np.array(indexies, dtype=np.int32), mode)
            else:
                write_array_to_file("lalal/indexies.bin", np.array(indexies[1:], dtype=np.int32), mode)

            write_array_to_file("lalal/vals.bin", np.array(vals, dtype=np.int32), mode)
            del indexies
            del vals

            print(j)
            j += len_of_pixels

        print("time: ", time.time() - start_time)


    def get_lines(self, chr : int, start : int, WINDOW_SIZE : int):
        window_size =  WINDOW_SIZE // 1_000
        return read_convert_to_hic_matrix(start // 1_000, window_size)

    def get_name(self):
        return "hi_c_limit"
        

        

    