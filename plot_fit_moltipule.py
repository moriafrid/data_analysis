import signal
from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt
from open_pickle import read_from_pickle
import os
from tqdm import tqdm
from add_figure import add_figure
import quantities as pq
from glob import glob
import pickle
from calculate_F_factor import calculate_F_factor
path_single_traces=glob('data/traces_img/2017_05_08_A_0006/*pA.p')
path=path_single_traces[0]
I=int(path[path.rfind('/')+1:path.rfind('pA')])
# path1='/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/data/short_pulse_mean.p'
path='/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/data/short_pulse/mean_short_pulse_with_parameters.p'

I=-50

SPINE_START = 60
resize_diam_by=1.5
CM=1#2/2
RM=5684*4#*2
RA=100

start_fit=0# 2960
end_fit=8900#3960
save_num='5'
do_calculate_F_factor=True
folder_ = "data/fit/"
try:os.mkdir(folder_)
except FileExistsError:pass
try:os.mkdir(folder_ +  'fits_together/')
except FileExistsError:pass

print('the injection current is',I,flush=True)
h.load_file("import3d.hoc")
h.load_file("nrngui.hoc")
h.load_file('stdlib.hoc')
h.load_file("stdgui.hoc")
# h.loadfile("stdrun.hoc")
if do_calculate_F_factor:
    V_head=0.14
    spine_neck_diam=0.164
    spine_neck_L=0.782



def SIGSEGV_signal_arises(signalNum, stack):
    print(f"{signalNum} : SIGSEGV arises")
    # Your code

signal.signal(signal.SIGSEGV, SIGSEGV_signal_arises)

#F_factor is the correction from the spine A(spine area)/area(without spine)
# F_factor = {} - take the F_factor for the specific segment
def change_model_pas(CM=1, RA = 250, RM = 20000.0, E_PAS = -70.0, F_factor = {}):
    #input the neuron property    h.dt = 0.1

    h.distance(0,0.5, sec=soma) # it isn't good beacause it change the synapse distance to the soma
    #h.distance(0, sec=soma)
    for sec in cell.all: ##check if we need to insert Ra,cm,g_pas,e_pas to the dendrit or just to the soma
        if isinstance(cell, Cell):
            if sec in cell.axon: continue   #@# cell.axon is not exist in hoc object
        sec.Ra = RA
        sec.cm = CM
        sec.g_pas = 1.0 / RM
        sec.e_pas = E_PAS
    for sec in cell.dend:
        for seg in sec: #count the number of segment and calclate g_factor and total dend distance,
            # how many segment have diffrent space larger then SPINE_START that decided
            if h.distance(seg) > SPINE_START:
                if do_calculate_F_factor:
                    F_factor=calculate_F_factor(cell,V_head,spine_neck_diam,spine_neck_L)
                else:
                    F_factor = 1.9#2.0 # F_factor[sec]
                seg.cm *= F_factor
                seg.g_pas *= F_factor



## e_pas is the equilibrium potential of the passive current
def plot_res(RM, RA, CM, save_name= "fit",print_full_graph=False,name='fit'):
    change_model_pas(CM=CM, RA=RA, RM=RM, E_PAS = E_PAS)
    Vvec = h.Vector() #cerat vector to record on it
    Tvec = h.Vector() #cerat vector to record on it
    Vvec.record(soma(0.5)._ref_v) #where to recprd
    Tvec.record(h._ref_t) #when it record
    h.cvode.store_events(Vvec)
    h.run()
    npTvec = np.array(Tvec)
    npVec = np.array(Vvec)
    plt.plot(npTvec[start_fit:end_fit], npVec[start_fit:end_fit],lw=2, alpha=0.5,linestyle ="-",color='green',label=name+' '+str([RM,RA,CM]) )#plot the recorded short_pulse
    exp_V = V#[int(180.0 / h.dt):int(800.0 / h.dt)]
    npVec = npVec#[int(180.0 / h.dt):int(800.0 / h.dt)]
    npVec = npVec[:len(exp_V)]
    error_1 = np.sqrt(np.sum(np.power(np.mean(exp_V[:2000]) - np.mean(npVec[:2000]), 2)))  # error from mean rest
    error_2 = np.sqrt(np.sum(np.power(exp_V[start_fit:end_fit] - npVec[start_fit:end_fit], 2))/(end_fit-start_fit))  #  error for the decay
    error_3 = np.sqrt(np.sum(np.power(np.mean(exp_V[4100:4900]) - np.mean(npVec[4100:4900]), 2)))  # error for maximal voltage
    error_tot = np.sqrt(np.sum(np.power(exp_V - npVec, 2))/len(exp_V)) # mean square error

    print('error_total=',round(error_tot,3))
    print('error_decay=', round(error_2,3))
    print('error_mean_max_voltage=', round(error_3,3))
    print('error_from_rest=', round(error_1,3))
    if print_full_graph:
        add_figure('fit'+str(I)+'pA part ['+str(start_fit)+':'+str(end_fit)+']',short_pulse[0].units,short_pulse[0].units)
        plt.plot(T, V, color = 'k') #plot short_pulse data
        plt.plot(npTvec[:len(npVec)], npVec, color = 'r', linestyle ="--") #plot the recorded short_pulse
        plt.plot(T[start_fit:end_fit], V[start_fit:end_fit], color = 'green')
        plt.suptitle('error from full graph='+str(round(error_tot,3))+' and eddor from decay='+str(round(error_2,3)))
        plt.legend(['data','NEURON_sim','decay_to_fitting'])
        plt.savefig(folder_+str(I)+'pA/'+save_name+"_full_graph.png")
    a=1

class Cell: pass
def mkcell(fname):
    #def to read ACS file
    h('objref cell, tobj')
    loader = h.Import3d_GUI(None)
    loader.box.unmap()
    loader.readfile(fname)
    c = Cell()
    loader.instantiate(c)
    return c
def instantiate_swc(filename):
    h('objref cell, tobj')
    h.load_file('allen_model.hoc')
    h.execute('cell = new allen_model()')
    h.load_file(filename)
    nl = h.Import3d_SWC_read()
    nl.quiet = 1
    nl.input(filename)
    i3d = h.Import3d_GUI(nl, 0)
    i3d.instantiate(h.cell)
    return h.cell

######################################################
# build the model
######################################################

fname = "05_08_A_01062017_Splice_shrink_FINISHED_LABEL_Bluecell_spinec91.ASC"
# cell=instantiate_swc('/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/try1.swc')
cell =mkcell(fname)
print (cell)
sp = h.PlotShape()
sp.show(0)  # show diameters

## delete all the axons
for sec in cell.axon:
   h.delete_section(sec=sec)

soma= cell.soma[0]

for sec in h.allsec():
    if sec == cell.soma[0]: continue
    sec.diam = sec.diam*resize_diam_by

#insert pas to all other section

for sec in tqdm(h.allsec()):
    sec.insert('pas') # insert passive property
    sec.nseg = int(sec.L/10)+1  #decide that the number of segment will be 21 with the same distances


# clamp = h.IClamp(cell.dend[82](0.996)) # insert clamp(constant potentientiol) at the soma's center
clamp = h.IClamp(soma(0.5)) # insert clamp(constant potentientiol) at the soma's center
clamp.amp = I/1000#-0.05 ## supopsed to be 0.05nA
if path in path_single_traces:
    start_inj=10500
    end_inj= 20400
    hz=0.1
    clamp.dur = (end_inj-start_inj)*hz
    clamp.delay = start_inj*hz
else:
    clamp.dur = 200
    clamp.delay = 296
######################################################
# load the data and see what we have
######################################################
short_pulse_dict = read_from_pickle(path)
short_pulse=short_pulse_dict['mean']
V = np.array(short_pulse[0])
short_pulse[1]=short_pulse[1].rescale('ms')
T = np.array(short_pulse[1])
T = T-T[0]
T=T
if path in path_single_traces:
    E_PAS = np.mean(V[:9000])
    end_fit=len(T)
else:
    E_PAS = np.mean(V[2945-500:2945])

h.tstop = (T[-1]-T[0])
h.v_init=E_PAS
h.dt = 0.1
h.steps_per_ms = h.dt

imp = h.Impedance(sec=soma)
imp.loc(soma(0.5))
imp.compute(0)
imp.input(0)
names=['paper par','free fit','CM<1','CM<1,RM>50','dend*2']
parameters=[[],[],[],[]]
dict={'paper par':[11368,100,1],'free fit':[12371,95.7,1.88],'CM=1':[12990.8,75.1,1],'RA=70':[13163.2,70,1.51],'RA=100':[12234.8,100,1.8],'RA=0':[15641,1e-9,1.2],'CM<1':[15595,1.0,0.97],'CM=1,RA=0':[15635.7,1e-9,1],'CM<1.2':[15579,1.5,1.2],'CM<1.2,RM>50':[13782.2,51,0.95]}
dict_after_resize={'1':[12392.7,95,1.86],'1.2':[15052.6,84.9,1.5],'1.5':[18057.4,93.2,1.13],'1.505':[14579.5,295.9,1.2],'1.51':[14639.3,298.7,1.71],'1.54':[15016.5,296,1.2],'1.55':[19199,59.6,1.01],'1.7':[20380,74.4,0.97],'1.8':[23035,87.8,0.88]}
if resize_diam_by!=1:
    name=str(resize_diam_by)
    add_figure("Fits [Rm,RA,Cm]", short_pulse[1].units, short_pulse[0].units)
    plt.plot(T[start_fit:end_fit], V[start_fit:end_fit], color='black', lw=4, alpha=0.2)
    RM, RA, CM=dict_after_resize[name]
    plot_res(RM, RA, CM,name='dend*'+name)
    imp.compute(0)
    print('the impadence is', imp.input(0))
    plt.legend(loc="lower left",prop={'size': 10})
    plt.savefig(folder_ +  'fits_together/dend*'+name+'.png')
    plt.savefig(folder_ +  'fits_together/dend*'+name+'.pdf')
    plt.close()
else:
    # for name in ['gregor par','free fit','RA=0','RA=70','CM=1,RA=0']:
    for name in dict.keys():
        add_figure("Fits [Rm,RA,Cm]", short_pulse[1].units, short_pulse[0].units)
        plt.plot(T[start_fit:end_fit], V[start_fit:end_fit], color='black', lw=4, alpha=0.2)
        RM, RA, CM=dict[name]
        plot_res(RM, RA, CM,name=name)
        imp.compute(0)
        print('the impadence is', imp.input(0))
        plt.legend(loc="lower left",prop={'size': 10})
        plt.savefig(folder_ +  'fits_together/'+name)
        plt.savefig(folder_ +  'fits_together/'+name+'.pdf')
        plt.close()


# add_figure("fit "+str(I)+"pA\nRM="+str(round(RM,1))+",RA="+str(round(RA,1))+",CM="+str(round(CM,2)),'mS','mV')












