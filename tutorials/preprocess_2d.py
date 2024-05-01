
from dnaloader.characteristics import Characteristic2dWithLimit
import os

if __name__ == '__main__':
    # You need to pass hic_path and bin_size to the function to pre-process a two-dimensional track
    current_path = os.getcwd()
    file_path_mcool = current_path + "/helper/4DNFI9FVHJZQ.mcool"
    track_2d = Characteristic2dWithLimit(
        folder_res="hic_result_1000",
        loader_type="hard",
        hic_path=file_path_mcool,
        bin_size=1_000)