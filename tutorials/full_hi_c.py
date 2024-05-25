from dnaloader.characteristics import CharacteristicFullHiC, CharacteristicHiCColer

if __name__ == '__main__':
    file_path_mcool = "/home/ojpochemy/SamplerBigWig/hi_c/4DNFIPO1DGLH.mcool"
    track_hic = CharacteristicFullHiC("/home/ojpochemy/dnaloader/hic_r") 
    print(track_hic.get_lines(chr=0, start_0=10000, start_1=9000, window_size=10000))
    # track_hic.preprocess()

    track_2d_coller = CharacteristicHiCColer(
        hic_path=file_path_mcool,
        bin_size=1_000
    )

    # print(track_2d_coller.get_lines(chr=0, start=1000, WINDOW_SIZE=10000))

    print(track_2d_coller.get_line_d(0, 10000, 9000, 10000))

    