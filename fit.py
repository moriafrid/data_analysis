'3'
import signal
from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt
from open_pickle import read_from_pickle
import os
from tqdm import tqdm
from add_figure import add_figure
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
shrinkage_factor=1#1.0/0.7
resize_diam_by=1

CM=1#2/2
RM=15000#5684*2#*2
RA=100

start_fit= 2000
end_fit=4900#3960

print('the injection current is',I,flush=True)
folder_='/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/'
h.load_file("import3d.hoc")
h.load_file("nrngui.hoc")
h.load_file('stdlib.hoc')
h.load_file("stdgui.hoc")
# h.loadfile("stdrun.hoc")
do_calculate_F_factor=True
spine_type="mouse_spine"
if spine_type=="groger_spine":
    V_head=0.14
    spine_neck_diam=0.164
    spine_neck_L=0.782
elif spine_type=="mouse_spine":
    #mouse
    from math import pi
    head_area=0.37
    r_head=np.sqrt(head_area/(4*pi))
    spine_neck_L=0.73
    # HEAD_DIAM=0.667
    spine_neck_diam=0.25 #0.164/07=0.23
    spine_density=1.08
    V_head=4/3*pi*r_head**3

def SIGSEGV_signal_arises(signalNum, stack):
    print(f"{signalNum} : SIGSEGV arises")
    # Your code
signal.signal(signal.SIGSEGV, SIGSEGV_signal_arises)
def get_inj(T,I,V):
    #found the begining,end and median of the injection
    I_abs = np.abs(I)
    inj_start = np.where(I_abs > I_abs.max() / 4.0)[0][0] - 1
    inj_end = np.where(I_abs > I_abs.max() / 4.0)[0][-1]
    inj = np.median(I[inj_start:inj_end])
    return inj, T[inj_start], T[inj_end]

# the function that looking for the best parameter fit
#F_factor is the correction from the spine A(spine area)/area(without spine)
# F_factor = {} - take the F_factor for the specific segment
def change_model_pas(CM=1, RA = 250, RM = 20000.0, E_PAS = -70.0, F_factor = {}):
    #input the neuron property    h.dt = 0.1

    h.distance(0,0.5, sec=soma) # it isn't good beacause it change the synapse distance to the soma
    #h.distance(0, sec=soma)
    for sec in h.allsec(): ##check if we need to insert Ra,cm,g_pas,e_pas to the dendrit or just to the soma
        sec.Ra = RA
        sec.cm = CM#*(1.0/0.7)
        sec.g_pas = (1.0 / RM)#*(1.0/0.7)
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
def plot_res(RM, RA, CM, save_name= "fit",print_full_graph=False):
    folder_="data/fit/"
    try: os.mkdir(folder_)
    except FileExistsError: pass
    if resize_diam_by!=1 and shrinkage_factor!=1:
        try: os.mkdir(folder_+'dend*'+str(round(resize_diam_by,2))+' &shrinkage by '+str(round(shrinkage_factor,2))+'_'+str(I)+'pA/')
        except FileExistsError: pass
        save_folder=folder_+'dend*'+str(round(resize_diam_by,2))+' &shrinkage by '+str(round(shrinkage_factor,2))+'_'+str(I)+'pA/'
    elif resize_diam_by!=1:
        try: os.mkdir(folder_+'dend*'+str(resize_diam_by)+'_'+str(I)+'pA/')
        except FileExistsError: pass
        save_folder=folder_+'dend*'+str(resize_diam_by)+'_'+str(I)+'pA/'
    elif shrinkage_factor!=1:
        try: os.mkdir(folder_+'F_shrinkage='+str(round(shrinkage_factor,2))+'_'+str(I)+'pA/')
        except FileExistsError: pass
        save_folder=folder_+'F_shrinkage='+str(round(shrinkage_factor,2))+'_'+str(I)+'pA/'
    else:
        try: os.mkdir(folder_+str(I)+'pA/')
        except FileExistsError: pass
        save_folder=folder_+str(I)+'pA/'

    # creat a clamp and record it for the chosen parameter
    ## save_name need to incloud the folder path
    change_model_pas(CM=CM, RA=RA, RM=RM, E_PAS = E_PAS)
    Vvec = h.Vector() #cerat vector to record on it
    Tvec = h.Vector() #cerat vector to record on it
    Vvec.record(soma(0.5)._ref_v) #where to recprd
    Tvec.record(h._ref_t) #when it record
    h.cvode.store_events(Vvec)

    h.run()
    npTvec = np.array(Tvec)
    npVec = np.array(Vvec)
    add_figure("fit "+str(I)+"pA\nRM="+str(round(RM,1))+",RA="+str(round(RA,1))+",CM="+str(round(CM,2)),'mS','mV')
    plt.plot(npTvec[start_fit:end_fit], npVec[start_fit:end_fit], color = 'r', linestyle ="--") #plot the recorded short_pulse
    plt.plot(T[start_fit:end_fit], V[start_fit:end_fit], color = 'green')
    plt.legend(['NEURON_sim','decay_to_fitting'])
    plt.savefig(save_folder+'/'+save_name+"_decay.png")
    plt.close()
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
        plt.savefig(save_folder+'/'+save_name+"_full_graph.png")
        plt.close()
    return save_folder




def efun(vals):
    #check the fitting
    # if the parameter incloud the fitting (not aqual to 1) check that the result is makes sense, if not return 1e6
    # if the result is make sense calculate the error between the record simulation and the initial data record
    ## *_IX is the parameter we play with them
    ## *_const is the basic parameters we return if the  result doesn't make sense
    if RM_IX != -1 :
        if vals.x[RM_IX] > 100000:
            return (1e6)
        RM = vals.x[RM_IX]
    else: RM = RM_const

    if CM_IX != -1:
        if vals.x[CM_IX] >2 :
            return (1e6)
        CM = vals.x[CM_IX]
    else:CM = CM_const

    if RA_IX != -1:
        if vals.x[RA_IX] > 300:
            return (1e6)
        RA = vals.x[RA_IX]
    else:RA = RA_const
    if (CM < 0.3 or RM < 2000 or RA <50):
        return 1e6
    # print('RA:',RA, '   CM:',CM, '   RM:',RM)

    change_model_pas(CM=CM, RA=RA, RM = RM, E_PAS = E_PAS)
    Vvec = h.Vector()
    Vvec.record(soma(0.5)._ref_v)

    h.run()
    npVec = np.array(Vvec)

    exp_V = V#[int(180.0/h.dt):int(800.0/h.dt)]
    npVec = npVec#[int(180.0/h.dt):int(800.0/h.dt)]
    npVec = npVec[:len(exp_V)]
    error_tot = np.sqrt(np.sum(np.power(exp_V - npVec, 2)))#/len(exp_V)) # mean square error


    error_1 = np.sqrt(np.sum(np.power(np.mean(exp_V[:2000]) - np.mean(npVec[:2000]), 2)))  # error from mean rest
    error_2 = np.sqrt(np.sum(np.power(exp_V[start_fit:end_fit] - npVec[start_fit:end_fit], 2))) #/(end_fit-start_fit)  #  error for the decay
    error_3 = np.sqrt(np.sum(np.power(np.mean(exp_V[4100:4900]) - np.mean(npVec[4100:4900]), 2)))  # error for maximal voltage

    return error_2 + (end_fit-start_fit)*error_3  #@# ask yoni if the calculation is right

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
#########################################
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

#insert pas to all other section
for sec in tqdm(h.allsec()):
    sec.insert('pas') # insert passive property
    sec.nseg = int(sec.L/10)+1  #decide that the number of segment will be 21 with the same distances
for sec in h.allsec():
    sec.diam = sec.diam*resize_diam_by
    sec.L*=shrinkage_factor



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
    # E_PAS = np.mean(V[:2000])

h.tstop = (T[-1]-T[0])
h.v_init=E_PAS
h.dt = 0.1
h.steps_per_ms = h.dt

CM_IX = 2
RM_IX=0
RA_IX = 1

RM_const = 60000.0
RA_const = 100
CM_const = 1.0

print("free params:")

h.attr_praxis(1e-9,1000,0)
opt_vals = h.Vector(3)#3
opt_vals.x[RM_IX] =RM
opt_vals.x[RA_IX] = RA
opt_vals.x[CM_IX] = CM
change_model_pas(CM=CM, RA=RA, RM=RM, E_PAS=E_PAS)

imp = h.Impedance(sec=soma)
imp.loc(soma(0.5))
imp.compute(0)
imp.input(0)

plot_res(CM=CM, RM=RM, RA=RA, save_name="before")
# plot_res(CM=CM, RM=RM, RA=RA, save_name="before")

#plot_res(CM=0.97, RM=11853.6, RA=99.6, save_name="test")

print('the initial impadence is', imp.input(0))

# allway run the fitting 3 time to avoid stack in local minima
for i in range(3):
    RMSD = h.fit_praxis(efun,opt_vals)   #@# take too much time if the fitting isn't found
    print('RMSD:',RMSD)
    RM = opt_vals.x[RM_IX]
    RA = opt_vals.x[RA_IX]
    CM = opt_vals.x[CM_IX]

    print("RMSD", RMSD,", RM",  RM, ", RA",  RA, ", CM",  CM)
    if i==2:
        plot_res(CM=CM, RM=RM, RA=RA, save_name="_fit_after_" + str(i + 1), print_full_graph=True)
    else:
        save_folder=plot_res(CM=CM, RM=RM, RA=RA, save_name="_fit_after_" + str(i + 1))

    # imp.compute(0)
    # print('the impadence is',imp.input(0))
pickle.dump({
        "RM": RM,
        "RA": RA,
        "CM": CM
    }, open(save_folder+'/' + "final_result_dend*"+str(resize_diam_by)+".p", "wb"))
#






