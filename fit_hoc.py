import signal
from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt
from open_pickle import read_from_pickle
import os

SPINE_START = 60 # mabe better to change to ~0._
hz=10000

path='/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/data/short_pulse_mean.p'
folder_='/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/'
h.load_file("import3d.hoc")
h.load_file("nrngui.hoc")
h.load_file('stdlib.hoc')
h.load_file("stdgui.hoc")



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
    for sec in h.allsec():
        sec.Ra = RA
        sec.cm = CM
        sec.g_pas = 1.0 / RM
        sec.e_pas = E_PAS
    for sec in h.dend:
        for seg in sec: #count the number of segment and calclate g_factor and total dend distance,
            # how many segment have diffrent space larger then SPINE_START that decided
            if h.distance(seg) > SPINE_START:
                F_factor = 2 # F_factor[sec]
                seg.cm *= F_factor
                seg.g_pas *= F_factor



## e_pas is the equilibrium potential of the passive current
def plot_res(RM, RA, CM, save_name= "fit"):
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
    Tvec.record(h._ref_t) #where to recprd

    h.run()
    npTvec = np.array(Tvec)
    npVec = np.array(Vvec)
    plt.plot(T, V, color = 'k') #plot short_pulse data
    plt.plot(npTvec, npVec, color = 'r', linestyle ="--") #plot the recorded short_pulse
    plt.title("fit\nRM="+str(round(RM,1))+",RA="+str(round(RA,1))+",CM="+str(round(CM,2)))
    plt.xlabel("time (ms)")
    plt.ylabel("V (mvdt)")
    plt.savefig(folder_+save_name+"_hoc.png")
    #plt.savefig(folder_+save_name+"_long.pdf")
    plt.close()
    exp_V = V#[int(180.0 / h.dt):int(800.0 / h.dt)]
    npVec = npVec#[int(180.0 / h.dt):int(800.0 / h.dt)]
    npVec = npVec[:len(exp_V)]
    error_2 = np.sqrt(np.sum(np.power(exp_V - npVec, 2)))  # mean square error
    print('error=', error_2)


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
        if vals.x[CM_IX] > 5 :
            return (1e6)
        CM = vals.x[CM_IX]
    else:CM = CM_const

    if RA_IX != -1:
        if vals.x[RA_IX] > 350:
            return (1e6)
        RA = vals.x[RA_IX]
    else:RA = RA_const

    #if (CM < 0.3 or RM < 10000 or RA < 100):
    #    return 1e6
    if (CM < 0.3 or RM < 2000 or RA < 1):
        return 1e6

    change_model_pas(CM=CM, RA=RA, RM = RM, E_PAS = E_PAS)
    Vvec = h.Vector()
    Vvec.record(soma(0.5)._ref_v)

    h.run()
    npVec = np.array(Vvec)
    exp_V = V#[int(180.0/h.dt):int(800.0/h.dt)]
    npVec = npVec#[int(180.0/h.dt):int(800.0/h.dt)]
    npVec = npVec[:len(exp_V)]
    error_2 = np.sqrt(np.sum(np.power(exp_V - npVec, 2))) # mean square error
    return error_2



######################################################
# build the model
######################################################


sp = h.PlotShape()
sp.show(0)  # show diameters


#post_morph_dir = "171101HuSHS2C1IN0toIN1__postsynaptic_reconstruction_with_putative_synapses.ASC"

model_file = "allen_model"
model_file = "morphology_correctedZ"
model_path = ""
h.load_file(model_file + ".hoc")
sp = h.PlotShape()
sp.show(0)  # show diameters
# delete all the axons
for sec in h.axon:
   h.delete_section(sec=sec)
soma = h.soma

# insert pas to all other section
for sec in h.allsec():
    sec.insert('pas')
    sec.nseg = int(sec.L/10)+1  #decide that the number of segment will be 21 with the same distances

clamp = h.IClamp(soma(0.5)) # insert clamp(constant potentientiol) at the soma's center
clamp.amp = -0.05
clamp.dur = 200
clamp.delay = 190
######################################################
# load the data and see what we have
######################################################
short_pulse, hz, rest = read_from_pickle(path, hz=True ,rest=True)
V = short_pulse[0]
T = short_pulse[1]
E_PAS = rest
V+=E_PAS

h.tstop = T[-1]
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
CM=2.7
RM=14893.5
RA=26.8

#CM=1
#RM=60000
#RA=200
opt_vals.x[RM_IX] =RM
opt_vals.x[RA_IX] = RA
opt_vals.x[CM_IX] = CM

imp = h.Impedance(sec=soma)
imp.loc(soma(0.5))
imp.compute(0)
imp.input(0)

plot_res(CM=CM, RM=RM, RA=RA, save_name="before")

# allway run the fitting 3 time to avoid stack in local minima
for i in range(3):
    RMSD = h.fit_praxis(efun,opt_vals)

    RM = opt_vals.x[RM_IX]
    RA = opt_vals.x[RA_IX]
    CM = opt_vals.x[CM_IX]

    print("RMSD", RMSD,", RM",  RM, ", RA",  RA, ", CM",  CM)
    plot_res(CM=CM, RM=RM, RA = RA, save_name="_fit_after_"+str(i+1))
    print("RMSD", RMSD,", RM",  RM, ", RA",  RA, ", CM",  CM)
    plot_res(CM=CM, RM=RM, RA = RA, save_name="_fit_after_"+str(i+1))
    imp.compute(0)
    print('the impadence is',imp.input(0))

soma_ref=h.SectionRef(sec=h.soma)
for i in range(soma_ref.nchild()):
    print(soma_ref.child[i](0).diam)
### othe function:

# h.dt = M1[:,0][1]-M1[:,0][0]
# h.steps_per_ms = 1.0/h.dt
# h.tstop = inj1_end + 100
# V_INIT = E_PAS
# h.v_init = E_PAS



# change_model_pas()

# stim = hcell.IClamp(soma(0.5))
# stim.delay = 50
# stim.dur = 2
# stim.amp = 0.075
#
# soma_v = h.Vector()
# soma_v.record(soma(0.5)._ref_v)
# I = h.Vector()
# I.record(stim._ref_i)
# t = h.Vector()
# t.record(h._ref_t)
# h.v_init = -70
# h.tstop = 200
# h.run()
#
# import numpy as np
# # np.savetxt("short_pll
# import matplotlib.pyplot as plt
# plt.plot(t,soma_v)
# plt.figure()
# plt.plot(t,I)
# plt.show()
# a=1

