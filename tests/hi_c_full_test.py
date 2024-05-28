import random
import numpy as np
import os
from dnaloader.characteristics import CharacteristicFullHiC, CharacteristicHiCColer

if __name__ == '__main__':
    # TODO: rewrite in pytest
    current_path = os.getcwd()
    file_path = "/home/ojpochemy/SamplerBigWig/hi_c/4DNFIPO1DGLH.mcool"

    chr_full = CharacteristicFullHiC("/home/ojpochemy/dnaloader/hic_test.h5")
    chr_with_coller = CharacteristicHiCColer(
        file_path, 1_000)

    for i in range(0, 1000):
        random_chr = 0
        chr_size = chr_full.get_chr_lenght(random_chr)
        random_number_1 = random.randint(0, chr_size)
        random_number_2 = random.randint(0, chr_size)
        random_window = random.randint(0, min(1_000_000, chr_size - max(random_number_1, random_number_2)))

        my_answer = chr_full.get_lines(
            random_chr,
            random_number_1,
            random_number_2,
            random_window)

        right_answer = chr_with_coller.get_line_d(random_chr, random_number_1, random_number_2, random_window)
        if (len(my_answer) == len(right_answer) and 
            len(my_answer[0]) == len(right_answer[0])):
            print("Is a square")
            if not np.array_equal(np.array(my_answer), np.array(right_answer)):
                raise ValueError(
                    "The two arrays are not the same! The program has been stopped.",
                    random_number_1,
                    " ",
                    random_number_2,
                    " ",
                    random_window)
        else:
            print("Not a square")
