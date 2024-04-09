from .Characteristics import Сharacteristic
import cooler

class СharacteristicHiCColer(Сharacteristic):
    def __init__(self, hic_path : str):
        super().__init__(
            hic_path,
        )
        self.c_matrix = cooler.Cooler(f"{hic_path}::/resolutions/5000")

    def get_lines(self, chr : int, start : int, WINDOW_SIZE : int):
        return self.c_matrix.matrix(balance=True).fetch((self.get_chr_name(chr + 1), start, start + WINDOW_SIZE))

    def get_name(self):
        return "hi_c"

    