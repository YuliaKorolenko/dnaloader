
from dnaloader.characteristics import Characteristic2dWithLimit
import os

if __name__ == '__main__':
    # You need to pass hic_path and bin_size to the function to pre-process a two-dimensional track
    current_path = os.getcwd()
    # file_path_mcool - the path to the track to be preprocessed
    # res - the path to the folder where the preprocessing result will be located
    file_path_mcool = current_path + "/helper/4DNFI9FVHJZQ.mcool"
    res = "hic_result_1000"

    track_2d = Characteristic2dWithLimit(
        folder_res=res,
        loader_type="hard",
        hic_path=file_path_mcool,
        bin_size=1_000)