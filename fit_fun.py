import signal
from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt
from open_pickle import read_from_pickle
import os
from tqdm import tqdm
from add_figure import add_figure
import quantities as pq

SPINE_START = 60
resize_diam_by=1
hz=10000

start_fit=2960
end_fit=3960#3960 #3105

path1='/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/data/short_pulse_mean.p'
path='/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/data/short_pulse/mean_short_pulse.p'
folder_='/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/'
h.load_file("import3d.hoc")
h.load_file("nrngui.hoc")
h.load_file('stdlib.hoc')
h.load_file("stdgui.hoc")
# h.loadfile("stdrun.hoc")




# useto open the ASc file, for now isn't working
class MyCell:
    def __init__(self):
        morph_reader = h.Import3d_Neurolucida3()
        morph_reader.input('/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/05_08_A_01062017_Splice_shrink_FINISHED_LABEL_Bluecell_spinec91.ASC')
        i3d = h.Import3d_GUI(morph_reader, 0)
        i3d.instantiate(self)
#m = MyCell()


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
                F_factor = 2.0 # F_factor[sec]
                seg.cm *= F_factor
                seg.g_pas *= F_factor



## e_pas is the equilibrium potential of the passive current
def plot_res(RM, RA, CM, save_name= "fit",print_full_graph=False):
    folder_="data/fit/"
    try: os.mkdir(folder_)
    except FileExistsError: pass
    # function to plot the data and create 2 pulses with diffrent e_pas
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
    if print_full_graph:
        add_figure('fit_part ['+str(start_fit)+':'+str(end_fit)+']','mS','mV')
        plt.plot(T*1000, V, color = 'k') #plot short_pulse data
        plt.plot(npTvec, npVec, color = 'r', linestyle ="--") #plot the recorded short_pulse
        plt.plot(T[start_fit:end_fit]*1000, V[start_fit:end_fit], color = 'green')
        plt.legend(['data','NEURON_sim','decay_to_fitting'])
        plt.savefig(folder_+save_name+"_full_graph.png")
        plt.close()
    add_figure("fit\nRM="+str(round(RM,1))+",RA="+str(round(RA,1))+",CM="+str(round(CM,2)),'mS','mV')
    plt.plot(npTvec[start_fit:end_fit], npVec[start_fit:end_fit], color = 'r', linestyle ="--") #plot the recorded short_pulse
    plt.plot(T[start_fit:end_fit]*1000, V[start_fit:end_fit], color = 'green')
    plt.legend(['NEURON_sim','decay_to_fitting'])
    plt.savefig(folder_+save_name+"_decay.png")
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

    a=1




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
        if vals.x[CM_IX] >5 :
            return (1e6)
        CM = vals.x[CM_IX]
    else:CM = CM_const

    if RA_IX != -1:
        if vals.x[RA_IX] > 350:
            return (1e6)
        RA = vals.x[RA_IX]
    else:RA = RA_const

    if (CM < 0.3 or RM < 2000 or RA <50):
        return 1e6

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
    error_2 = np.sqrt(np.sum(np.power(exp_V[start_fit:end_fit] - npVec[start_fit:end_fit], 2))/(end_fit-start_fit))  #  error for the decay
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

######################################################
# build the model
######################################################

fname = "05_08_A_01062017_Splice_shrink_FINISHED_LABEL_Bluecell_spinec91.ASC"
# cell=instantiate_swc('/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/try1.swc')
def fit(cell_file,trace):
    cell =mkcell(fname)
    print (cell)
    sp = h.PlotShape()
    sp.show(0)  # show diameters

    ## delete all the axons
    for sec in cell.axon:
       h.delete_section(sec=sec)

    soma= cell.soma[0]
    for sec in soma.children():
        if isinstance(cell, Cell):
            if sec in cell.axon: continue
        sec.diam*=resize_diam_by

    #insert pas to all other section

    for sec in tqdm(cell.all):
        if isinstance(cell, Cell):
            if sec in cell.axon: continue   #@# cell.axon is not exist in hoc object
        sec.insert('pas') # insert passive property
        sec.nseg = int(sec.L/10)+1  #decide that the number of segment will be 21 with the same distances


    clamp = h.IClamp(soma(0.5)) # insert clamp(constant potentientiol) at the soma's center
    clamp.amp = -0.05 ## supopsed to be 0.05nA
    clamp.dur = 200
    clamp.delay = 296
    ######################################################
    # load the data and see what we have
    ######################################################
    short_pulse = read_from_pickle(path)
    V = np.array(short_pulse[1])
    T = np.array(short_pulse[0])
    T = T-T[0]
    E_PAS = np.mean(V[:2000])
    #V+=E_PAS

    h.tstop = (T[-1]-T[0])*1000
    h.v_init=E_PAS
    h.dt = 0.1
    h.steps_per_ms = h.dt
    CM_IX = 0
    RM_IX=1
    RA_IX = 2

    RM_const = 60000.0
    RA_const = 250.0
    CM_const = 1.0

    print("free params:")

    h.attr_praxis(1e-9,1000,0)
    opt_vals = h.Vector(3)
    CM=2/2
    RM=5684*2
    RA=100
    opt_vals.x[RM_IX] =RM
    opt_vals.x[RA_IX] = RA
    opt_vals.x[CM_IX] = CM
    change_model_pas(CM=CM, RA=RA, RM=RM, E_PAS=E_PAS)

    imp = h.Impedance(sec=soma)
    imp.loc(soma(0.5))
    imp.compute(0)
    imp.input(0)

    plot_res(CM=CM, RM=RM, RA=RA, save_name="before")
    #plot_res(CM=0.97, RM=11853.6, RA=99.6, save_name="test")

    print('the initial impadence is', imp.input(0))

    # allway run the fitting 3 time to avoid stack in local minima
    for i in range(3):
        RMSD = h.fit_praxis(efun,opt_vals)   #@# take too much time if the fitting isn't found

        RM = opt_vals.x[RM_IX]
        RA = opt_vals.x[RA_IX]
        CM = opt_vals.x[CM_IX]

        print("RMSD", RMSD,", RM",  RM, ", RA",  RA, ", CM",  CM)
        if i==2:
            plot_res(CM=CM, RM=RM, RA=RA, save_name="_fit_after_" + str(i + 1), print_full_graph=True)
        else:
            plot_res(CM=CM, RM=RM, RA=RA, save_name="_fit_after_" + str(i + 1))

        imp.compute(0)
        print('the impadence is',imp.input(0))

    soma_ref = h.SectionRef(sec=cell.soma[0])
    print("the soma's childrens diameter is:")
    for i in range(soma_ref.nchild()):
        print(soma_ref.child[i](0).diam,soma_ref.child[i](0))
    length = 0
    for dend in cell.dend:
        length += dend.L
    print("total dendritic length is ", length)


# track from the terminals to the soma
def track_one(cell,terminal):
    h.distance(0, 0.5, sec=cell.soma[0])
    sec = terminal
    dis = []
    diam = []
    while sec != cell.soma[0]:
        if isinstance(cell, Cell):
            if sec in cell.axon: continue
        dis.append(h.distance(sec.parentseg()))
        sec_ref = h.SectionRef(sec=sec)
        diam.append(sec.diam)
        sec = sec_ref.parent
    return np.array(dis), np.array(diam)


terminals = []
for sec in cell.dend:
    if len(sec.children()) == 0:
        terminals.append(sec)
plt.close()
add_figure('dendritic tree without constant branchs', 'distance from soma', 'diameter')
i = 0
for terminal in terminals:
    dis, diam = track_one(terminal)
    if round(np.mean(np.array(diam)),3)!=round(diam[-1],3):
        i += 1
        plt.plot(dis, diam, alpha=0.5)
plt.savefig('temp_groger2')
print(i,' graph is printed')




