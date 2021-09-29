from matplotlib import pyplot as plt
from open_pickle import read_from_pickle
import numpy as np
import scipy.fftpack
from scipy.signal import savgol_filter
import tkinter
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from add_figure import add_figure

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

def check_dynamics(short_pulse,save_folder):
    add_figure('smooth short_pulse','t[ms]','V[mV]')
    plt.plot(short_pulse,'.')
    yhat = savgol_filter(short_pulse, 25,2)  # window size 51, polynomial order 3
    plt.plot(yhat,'.')
    plt.savefig(save_folder+'/smooth_curve')
    plt.show()
    add_figure('2 part of short_pulse on each other','t[ms]','V[mV]')
    min_index=np.argmin(yhat)
    max_index=np.argmax(yhat)
    part1=yhat[max_index:min_index]-yhat[max_index]
    part2=-yhat[min_index:]+yhat[min_index]
    plt.plot(part1)
    plt.plot(part2)
    plt.savefig(save_folder+'/check_dynamics')
    plt.show()
    a=1



if __name__=='__main__':
    folder_short_pulse = '/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis/data/short_pulse/'
    a=read_from_pickle('/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis/data/short_pulse_correct.p')
    mean_short_pulse=np.mean(a,axis=0)
    check_dynamics(mean_short_pulse,folder_short_pulse)
    a=1
