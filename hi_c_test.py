import random
import numpy as np
from characteristics import CharacteristicHiCColer, CharacteristicHiCWithLimit

if __name__ == '__main__':
    # TODO: переписать в нормальный тест
    for i in range(0, 1000):
        random_number = random.randint(0, 220_000_000)
        random_window = random.randint(0, 1_000_000)

        chr_with_limit = CharacteristicHiCWithLimit("/home/ojpochemy/SamplerBigWig/hi_c/4DNFI9FVHJZQ.mcool", "/home/ojpochemy/dnaloader/hic_result")
        my_answer = chr_with_limit.get_lines(0, random_number, random_window)
        print("my answer: ", my_answer)

        chr_with_coller = CharacteristicHiCColer("/home/ojpochemy/SamplerBigWig/hi_c/4DNFI9FVHJZQ.mcool")
        right_answer = chr_with_coller.get_lines(0, random_number, random_window)

        print("right answer: ", right_answer)
        
        print(random_number)
        if not np.array_equal(np.array(my_answer), np.array(right_answer)):
            raise ValueError("Два массива не одинаковые! Программа остановлена. ", random_number, " ", random_window)
        