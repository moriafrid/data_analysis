from open_pickle import read_from_pickle
import numpy as np
import scipy.fftpack
import matplotlib.pyplot as plt
from add_figure import add_figure
import os

def correct_frequency(pheno):
    hz=10000
    N=len(pheno)
    x= np.arange(0, len(pheno), 1)*(1000.0/hz)
    y = pheno
    w = scipy.fftpack.rfft(y)
    f = scipy.fftpack.rfftfreq(N, x[1] - x[0])
    spectrum = w ** 2

    cutoff_idx = spectrum < (spectrum.max() / 2)
    w2 = w.copy()
    w2[cutoff_idx] = 0
    y2 = scipy.fftpack.irfft(w2)
    return y2

def check_dynamics(short_pulse,x_short_pulse,folder):
    try: os.mkdir(folder+'check_dynamic/')
    except FileExistsError: pass
    save_folder=folder+'check_dynamic/'

    add_figure('short_pulse',x_short_pulse.units,short_pulse.units)
    plt.plot(x_short_pulse,short_pulse,'.')
    #plt.plot(x_short_pulse[::2],short_pulse[::2],'.')
    plt.savefig(save_folder+'/short_pulse')
    plt.show()
    add_figure('2 part of short_pulse on each other',x_short_pulse.units/1000,short_pulse.units)
    min_index=np.argmin(short_pulse)
    max_index=np.argmax(short_pulse)
    part1=short_pulse[max_index:min_index]-short_pulse[max_index]
    part2=-short_pulse[min_index:]+short_pulse[min_index]
    x_part1=x_short_pulse[max_index:min_index]-x_short_pulse[max_index]
    x_part2=x_short_pulse[min_index:]-x_short_pulse[min_index]
    plt.plot(x_part1*1000, part1)
    plt.plot(x_part2[:len(x_part1)]*1000,part2[:len(x_part1)])
    plt.legend(['beginigng','flip end'])
    plt.savefig(save_folder+'/check_dynamics')
    plt.show()
    a=1



if __name__=='__main__':
    folder_short_pulse = '/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/data/check_dynamic/'
    mean_short_pulse,x_mean_short_pulse=read_from_pickle('/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/data/short_pulse/mean_short_pulse.p')
    check_dynamics(mean_short_pulse,x_mean_short_pulse,folder_short_pulse)
    a=1
