from neo import io
from matplotlib import pyplot as plt
import numpy as np
from tqdm import tqdm
from add_figure import add_figure
import pickle
from scipy.optimize import curve_fit
import quantities as pq


def linear(x, m):
    return m*x

def sepereat_by_current(t,T,I,f):
    add_figure('I_V curve_together', 'points', t.units)
    for v in np.array(t):
        plt.plot(v)
    plt.savefig('data/traces_img/' + f + '/' + 'I_V_curve_together')
    plt.close()
    maxi=[]
    for i,v in enumerate(np.array(t)):
        add_figure('Trace with current' + str(I[i]),'points',t.units)
        plt.plot(v,alpha=0.5)
        plt.plot(np.arange(len(v))[11600:20400], v[11600:20400])
        max_temp=np.mean(v[11600:20400])
        plt.plot([11600,20400],[max_temp,max_temp])
        rest=np.mean(t[I.index(0)][0:10000])
        plt.plot(np.arange(len(v))[0:10000], v[0:10000])
        rest_temp = np.mean(v[0:10000])
        plt.plot([0,10000], [rest_temp, rest_temp])
        plt.legend(['full trace','max to calculate','max'])
        plt.savefig('data/traces_img/'+f+'/'+str(I[i])+'pA.png')
        with open('data/traces_img/'+f+'/'+str(I[i])+'pA.p', 'wb') as file:
            pickle.dump([v*t.units,T[i]], file)
        plt.close()
        maxi.append(max_temp)
        rest=np.mean(t[I.index(0)][11600:20400])
    rest1=np.mean(t[I.index(0)])
    rest2=np.mean(t[I.index(0)][11600:20400])
    return (maxi-np.array(rest2))

def I_V_curve(maxi,I,save_file):
    add_figure('I-V Curve\nfit to V=I*Rinput','Current I[pA]','Voltage V[mV]')
    plt.plot(I,maxi,'.',label='max volt to diffrent current inject')
    popt, pcov = curve_fit(linear, np.array(I), np.array(maxi))
    plt.plot(I, linear(np.array(I), *popt),label='fit=I*'+str(round(popt[0]*1e-12/1e-3*1e12,3))+'pohm')
    plt.plot(I[-1],maxi[-1],'*',label=str(I[-1])+'pA')
    I_50=[0,I[-1]]
    maxi_50=[0,maxi[-1]]
    popt1, pcov1 = curve_fit(linear, np.array(I_50), np.array(maxi_50))
    plt.plot(I, linear(np.array(I), *popt1),label='fit=I*'+str(round(popt1[0]*1e-12/1e-3*1e12,3))+'pohm')
    plt.legend()
    # plt.legend(['max volt to diffrent current inject','fit=I*'+str(round(popt[0]*1e-12/1e-3*1e12,3))+'pohm',str(I[-1])+'pA'])
    plt.savefig( save_file+ 'I_V_curve_fit')
    print('The input resistance from I_V cureve is ',round(popt[0]*1e-12/1e-3*1e12,3),'pOhm')
    print('The input resistance from I=-50pA is ', round(popt1[0] * 1e-12 / 1e-3 * 1e12, 3), 'pOhm')
    return popt1[0]*10e-12/10e-3*pq.ohm
def find_maxi(V,f):
    plt.plot(V)
    plt.plot(np.arange(np.argmax(V)+750,np.argmin(V)-10),V[np.argmax(V)+750:np.argmin(V)-10])
    plt.plot(np.arange(np.argmax(V)+750,np.argmin(V)-10),np.mean(V[np.argmax(V)+750:np.argmin(V)-10])*np.ones((np.argmin(V)-10)-(np.argmax(V)+750)))
    plt.legend(['full trace','max to calculate','max'])
    plt.savefig('data/traces_img/' + f + '/-50pA.png')
    return np.mean(V[np.argmax(V)+550:np.argmin(V)-10])
if __name__=='__main__':
    I=[-200,-160,-120,-80,-40,-0,40,80,120,160]
    #    creat_data
    f='2017_05_08_A_0006'
    r=io.AxonIO(f+'.abf')
    bl = r.read_block(lazy=False)
    t_total,T_total,max_t=[],[],[]
    for i in tqdm(range(len(bl.segments))):
        segment =bl.segments[i]
        hz = segment.analogsignals[0].sampling_rate
        t_start=segment.analogsignals[0].t_start
        t_stop=segment.analogsignals[0].t_stop
        t = np.array(segment.analogsignals[0])
        V_unit=segment.analogsignals[0].units
        #t = np.array(t).flatten()
        T = np.linspace(t_start,t_stop, len(t))
        t_total.append(t[:,0])
        T_total.append(T)
    maxi=sepereat_by_current(t_total*V_unit,T_total,I,f)
    from open_pickle import read_from_pickle
    mean_short_pulse= read_from_pickle('/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/data/short_pulse/mean_short_pulse_with_parameters.p')
    t_50pA,T_50pA=mean_short_pulse['mean']
    t_50pA=np.array(t_50pA)-mean_short_pulse['E_pas']
    maxi=np.append(maxi,find_maxi(t_50pA,f))
    I.append(-50)
    # I_V_curve(t_total*V_unit, T_total,maxi,I,f)
    save_folder_IV_curve =  'data/traces_img/' + f + '/'

    I_V_curve(maxi, I, save_folder_IV_curve)
    plt.show()
    a=1

