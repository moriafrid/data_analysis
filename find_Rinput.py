#bring the data
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
from open_pickle import read_from_pickle
from add_figure import add_figure
from glob import glob
import numpy as np


def linear(x, m, c):
    return m*x+c

def find_Rinput(folder_):
    Rinput_list=[]
    for file in glob(folder_+'*.p'):
        t = read_from_pickle(file)
        T = np.linspace(-220, 220, len(t))
        t=t - np.mean(t)
        peaks, properties = find_peaks(abs(t), threshold=3, distance=20000)

        add_figure('I-V curve  - the data+fit','I [pA]','V [mV]')
        #plt.figure('seaborn')
        plt.plot(T, t )
        plt.plot(T[peaks],t[peaks],  "x")
        I=np.array([-200,-160,-120,-80,-40,80,120,160])
        #v=np.append(t[peaks],0)
        V = t[peaks]
        plt.plot(I,V,  "x")
        popt, pcov = curve_fit(linear, I, V)
        plt.plot(I, linear(I, *popt))
        #sns.regplot(peaks, T[peaks], "x")
        plt.savefig(file[:-2]+'data_and_fit')

        add_figure('fit for I-V curve','I [pA]','V [mV]')
        V = t[peaks]
        I2 = I[:,np.newaxis]
        #Rin, c,_, _ = np.linalg.lstsq(I2, V)
        plt.plot(I,V,  "x")
        plt.plot(0,0,'rx')
        plt.plot(I, linear(I, *popt),'b')
        plt.plot(I, popt[0]*I,'g')
        #plt.plot(I, Rin*I)
        Rinput_list.append(popt[0])
        plt.savefig(file[:-2]+'fit')

    print('the fit for the slope is '+str(popt[0])+' mV/pA',flush=True)
    print('the bias before constraint of [0,0] is: '+str(popt[1]),flush=True)
    print('the covarience is '+str(pcov),flush=True)
    plt.show()
    return np.mean(Rinput_list),Rinput_list

if __name__=='__main__':
    folder_ = '/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis/data/traces_img/2017_05_08_A_0006/'
    R,Rinput_list=find_Rinput(folder_)
    a=1