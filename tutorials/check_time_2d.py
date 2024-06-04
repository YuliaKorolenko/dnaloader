
from dnaloader.characteristics import Characteristic2dWithLimit, CharacteristicCooler, Characteristic
import os
import random
import time


def get_time(track: Characteristic, chr: int, start: int, window: int):
    start_time = time.time()
    track.get_lines(chr, start, window)
    return time.time() - start_time


if __name__ == '__main__':
    # You need to pass hic_path and bin_size to the function to pre-process a
    # two-dimensional track
    current_path = os.getcwd()
    # file_path_mcool - the path to the cool/mcool file
    # res - the path to the cool/mcool file

    file_path_mcool = current_path + "/helper/4DNFI9FVHJZQ.mcool"
    res = "hic_result_1000"

    track_2d = Characteristic2dWithLimit(
        folder_res=res,
        loader_type="hard")

    track_2d_coller = CharacteristicCooler(
        hic_path=file_path_mcool,
        bin_size=1_000
    )

    # Choose a random chromosome, starting position and window size
    random_chr = random.randint(0, 23)
    chr_size = track_2d.get_chr_lenght(random_chr)
    random_number = random.randint(0, chr_size)
    random_window = random.randint(0, min(1_000_000, chr_size - random_number))

    # Measure the time
    print(
        "2d with limit: ",
        get_time(
            track_2d,
            random_chr,
            random_number,
            random_window))
    print(
        "2d coller: ",
        get_time(
            track_2d_coller,
            random_chr,
            random_number,
            random_window))
