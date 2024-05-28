import random
import numpy as np
import os
from dnaloader.characteristics import CharacteristicHiCColer, Characteristic2dWithLimit

if __name__ == '__main__':
    # TODO: rewrite in pytest
    current_path = os.getcwd()
    file_path = current_path + "/helper/4DNFI9FVHJZQ.mcool"

    chr_with_limit = Characteristic2dWithLimit(
        "hic_result_1000", "hard")
    chr_with_coller = CharacteristicHiCColer(
        file_path, 1_000)

    for i in range(0, 10):
        random_chr = random.randint(0, 23)
        chr_size = chr_with_limit.get_chr_lenght(random_chr)
        random_number = random.randint(0, chr_size)
        random_window = random.randint(
            0, min(1_000_000, chr_size - random_number))

        my_answer = chr_with_limit.get_lines(
            random_chr, random_number, random_window)

        right_answer = chr_with_coller.get_lines(
            random_chr, random_number, random_window)
        if not np.array_equal(np.array(my_answer), np.array(right_answer)):
            raise ValueError(
                "The two arrays are not the same! The program has been stopped.",
                random_number,
                " ",
                random_window)
