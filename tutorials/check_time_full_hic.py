
from dnaloader.characteristics import CharacteristicFullHiC, CharacteristicHiCColer, Characteristic
import os
import random
import time
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as st

if __name__ == '__main__':
    current_path = os.getcwd()
    # file_path_mcool - the path to the cool/mcool file
    # res - the path to the cool/mcool file

    file_path_mcool = "/home/ojpochemy/SamplerBigWig/hi_c/4DNFIPO1DGLH.mcool"
    res = os.getcwd() + "/hic_2.h5"

    track_2d = CharacteristicFullHiC(res)

    track_2d_coller = CharacteristicHiCColer(
        hic_path=file_path_mcool,
        bin_size=1_000
    )

    data = {'My approach': [], 'Cooler': [], 'Window Size': []}
    for i in range(0, 4):

        cur_window_size = 1000 * (10 ** i)
        print("window size: ", cur_window_size)
        for i in range(0, 100):
            # Choose a random chromosome, starting position and window size
            random_chr = 0
            chr_size = track_2d.get_chr_lenght(random_chr) - cur_window_size
            random_1 = random.randint(0, chr_size)
            random_2 = random.randint(0, chr_size)

            start_time = time.time()
            a1 = track_2d.get_lines(
                random_chr, random_1, random_2, cur_window_size)
            data['My approach'].append(time.time() - start_time)

            start_time = time.time()
            a2 = track_2d_coller.get_lines(
                random_chr, random_1, random_2, cur_window_size)
            data['Cooler'].append(time.time() - start_time)

            data['Window Size'].append(cur_window_size)

    df = pd.DataFrame(data)

    sns.set(style="darkgrid")
    plt.figure(figsize=(9, 6))
    sns.lineplot(
        x='Window Size',
        y='My approach',
        data=df,
        label='My approach',
        errorbar=('ci'))
    sns.lineplot(
        x='Window Size',
        y='Cooler',
        data=df,
        label='Cooler',
        errorbar=('ci'))
    plt.xlabel('Размер окна')
    plt.ylabel('Время в секундах')
    plt.title('Время запроса для разного размера окна')
    plt.savefig('graph.png')

    for i in range(0, 4):
        cur_window_size = 1000 * (10 ** i)
        time1 = [data['My approach'][i] for i in range(
            len(data['My approach'])) if data['Window Size'][i] == cur_window_size]
        time2 = [data['Cooler'][i] for i in range(
            len(data['Cooler'])) if data['Window Size'][i] == cur_window_size]

        print(cur_window_size)
        # create 95% confidence interval for population mean weight
        res1 = st.t.interval(
            df=len(data) - 1,
            loc=np.mean(time1),
            scale=st.sem(time1),
            confidence=0.95)
        print("my res: ", res1, " mean: ", np.mean(time1))
        res2 = st.t.interval(
            df=len(data) - 1,
            loc=np.mean(time2),
            scale=st.sem(time2),
            confidence=0.95)
        print("coller res: ", res2, " mean: ", np.mean(time2))
