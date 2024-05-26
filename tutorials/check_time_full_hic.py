
from dnaloader.characteristics import CharacteristicFullHiC, CharacteristicHiCColer, Characteristic
import os
import random
import time

def get_time(track : Characteristic, chr : int, start_1 : int, start_2 : int, window : int, is_mine : bool):
    start_time = time.time()
    if (is_mine):
        track.get_lines(chr, start_1, start_2, window)
    else:
        track.get_line_d(chr, start_1, start_2, window)
    return time.time() - start_time


if __name__ == '__main__':
    # You need to pass hic_path and bin_size to the function to pre-process a two-dimensional track
    current_path = os.getcwd()
    # file_path_mcool - the path to the cool/mcool file
    # res - the path to the cool/mcool file
    
    file_path_mcool =  "/home/ojpochemy/SamplerBigWig/hi_c/4DNFIPO1DGLH.mcool"
    res = "hic_r"

    track_2d = CharacteristicFullHiC("/home/ojpochemy/dnaloader/hic_r4") 
    
    track_2d_coller = CharacteristicHiCColer(
        hic_path=file_path_mcool,
        bin_size=1_000
    )

    # Choose a random chromosome, starting position and window size
    random_chr = 0
    chr_size = 2000000
    random_number_1 = random.randint(0, chr_size)
    random_number_2 = random.randint(0, chr_size)
    random_window = random.randint(0, min(1_000_000, chr_size - random_number_1))

    # Measure the time
    print("2d full: ", get_time(track_2d, random_chr, random_number_1, random_number_2, random_window, True))
    print("2d coller: ", get_time(track_2d_coller, random_chr, random_number_1, random_number_2, random_window, False))

    


