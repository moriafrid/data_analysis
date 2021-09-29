from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt
import pickle
SPINE_START = 60
hz=10000
def get_inj(T,I,V):
    I_abs = np.abs(I)
    inj_start = np.where(I_abs > I_abs.max() / 4.0)[0][0] - 1
    inj_end = np.where(I_abs > I_abs.max() / 4.0)[0][-1]
    inj = np.median(I[inj_start:inj_end])
    return inj, T[inj_start], T[inj_end]
# the function that looking for the best parameter fit
#F_factor is the correction from the spine A(spine area)/area(without spine)
# F_factor = {} - take the F_factor for the specific segment
def change_model_pas(CM=1, RA = 250, RM = 20000.0, E_PAS = -70.0, F_factor = {}):
    h.distance(0, sec=soma)
    for sec in h.allsec():
        sec.Ra = RA
        sec.cm = CM
        sec.g_pas = 1.0 / RM
        sec.e_pas = E_PAS
    for sec in h.dend:
        for seg in sec:
            if h.distance(seg.x) > SPINE_START:
                F_factor = 1.5 # F_factor[sec]
                seg.cm *= F_factor
                seg.g_pas *= F_factor


# function to plot in the end
# creat a clamp and record it for the chosen parameter
def plot_res(RM, RA, CM, save_name= "fit"):
    change_model_pas(CM=CM, RA=RA, RM=RM, E_PAS = E_PAS1)
    clamp.amp = inj1
    clamp.dur = inj1_end - inj1_start
    clamp.delay = inj1_start
    Vvec1 = h.Vector()
    Vvec1.record(soma(0.5)._ref_v)
    h.tstop = M1[:,0][-1]
    h.v_init=E_PAS1
    h.run()
    npVec1 = np.array(Vvec1)
    plt.plot(M1[:,0], M1[:,2], color = 'k')
    time = np.arange(0,len(npVec1) * h.dt,h.dt)
    plt.plot(time, npVec1, color = 'r', linestyle ="--")
    plt.title("long pulse fit\nRM="+str(round(RM,1))+",RA="+str(round(RA,1))+",CM="+str(round(CM,2)))
    plt.savefig(save_name+"_long.png")
    plt.savefig(save_name+"_long.pdf")
    plt.close()

    change_model_pas(CM=CM, RA=RA, RM=RM, E_PAS = E_PAS2)
    clamp.amp = inj2
    clamp.dur = inj2_end - inj2_start
    clamp.delay = inj2_start
    Vvec2 = h.Vector()
    Vvec2.record(soma(0.5)._ref_v)
    h.tstop = M2[:,0][-1]
    h.v_init=E_PAS2
    h.run()
    npVec2 = np.array(Vvec2)

    plt.plot(M2[:,0], M2[:,2], color = 'k')
    time = np.arange(0,len(npVec2) * h.dt,h.dt)
    plt.plot(time, npVec2, color = 'r', linestyle ="--")
    plt.xlabel("time (ms)")
    plt.ylabel("V (mv)")
    plt.title("short pulse fit\nRM="+str(round(RM,1))+",RA="+str(round(RA,1))+",CM="+str(round(CM,2)))
    plt.savefig(save_name+"_short.png")
    plt.savefig(save_name+"_short.pdf")
    plt.close()
# if the parameter incloud the fitting (not aqual to 1) check that the result is makes sense, if not return 1e6
# if the result is make sense calcculate the error between the record simulation and the initial data record
def efun(vals):
    if RM_IX != -1 :
        if vals.x[RM_IX] > 100000:
            return (1e6)
        RM = vals.x[RM_IX]
    else: RM = RM_const

    if CM_IX != -1:
        if vals.x[CM_IX] > 2 :
            return (1e6)
        CM = vals.x[CM_IX]
    else:CM = CM_const
    if RA_IX != -1:
        if vals.x[RA_IX] > 350:
            return (1e6)
        RA = vals.x[RA_IX]
    else:RA = RA_const

    if (CM < 0.3 or RM < 5000 or RA < 100):
        return 1e6

    change_model_pas(CM=CM, RA=RA, RM = RM, E_PAS = E_PAS1)
    clamp.amp = inj1
    clamp.dur = inj1_end - inj1_start
    clamp.delay = inj1_start
    Vvec1 = h.Vector()
    Vvec1.record(soma(0.5)._ref_v)
    h.v_init=E_PAS1
    h.tstop = M1[:,0][-1]
    h.run()
    npVec1 = np.array(Vvec1)
    error_1 = abs(npVec1.max()- M1[:,2].max())

    change_model_pas(CM=CM, RA=RA, RM = RM, E_PAS = E_PAS2)
    clamp.amp = inj2
    clamp.dur = inj2_end - inj2_start
    clamp.delay = inj2_start
    Vvec2 = h.Vector()
    Vvec2.record(soma(0.5)._ref_v)
    exp_data = M2[:,2]

    h.v_init=E_PAS2
    h.tstop = M2[:,0][-1]
    h.run()
    npVec2 = np.array(Vvec2)
    fit_start_for_short_pulse = int((inj2_end+1)/h.dt)
    exp_data_fit_part = exp_data[fit_start_for_short_pulse:fit_start_for_short_pulse+int(100/h.dt)]
    sim_data_fit_part = npVec2[fit_start_for_short_pulse:fit_start_for_short_pulse+int(100/h.dt)]
    error_2 = np.sum(np.sqrt(np.power(exp_data_fit_part - sim_data_fit_part, 2))) # mean square error
    print("cur error =", error_1 + error_2," RM=" , RM, " RA=", RA, " CM=" , CM)
    return error_1 + error_2

######################################################
# build the model
######################################################

#post_morph_dir = "morph.swc"

##%% change the input file from swc to abf
##%% change to the correct morphology
h.load_file("import3d.hoc")
h.load_file("nrngui.hoc")
#h.load_file('stdlib.hoc')

post_morph_dir = "05_08_A_01062017_Splice_shrink_FINISHED_LABEL_Bluecell_spinec91.ASC"
post_morph_dir = "171101HuSHS2C1IN0toIN1__postsynaptic_reconstruction_with_putative_synapses.ASC"
morph_file = post_morph_dir
#model_file = "allen_model"
model_file = "morphology_correctedZ"
model_path = ""
h.load_file(model_file + ".hoc")

for sec in h.axon:
   h.delete_section(sec=sec)
soma = h.soma
for sec in h.allsec():
    sec.insert('pas')
    sec.nseg = int(sec.L/20)+1

clamp = h.IClamp(soma(0.5))

######################################################
# load the data and see what we have
######################################################
long_pulse = "long_pulse.txt"  ##%% I don't have long pulse
M1 = np.loadtxt(long_pulse)
M1[:,1] -= M1[:,1][0]
inj1, inj1_start, inj1_end = get_inj(M1[:,0], M1[:,1], M1[:,2])
E_PAS1 = np.mean(M1[:,2][:200]) # get e pas as the mean of the start of the trace

short_pulse = "sshort_pulsehort_pulse.txt"    ##%% change the file from text to pickle
M2 = np.loadtxt(short_pulse)
E_PAS2 = np.median(M2[:,2][:200]) # get e pas as the mean of the start of the trace
inj2, inj2_start, inj2_end = get_inj(M2[:,0], M2[:,1], M2[:,2])

### othe function:

# ploting the data
# fig, ax = plt.subplots(2,2)
# plt.subplots_adjust(hspace = 0.5, wspace = 0.5)
# ax[0][0].set_title("voltage long pulse")
# ax[0][0].plot(M1[:,0], M1[:,2], color="k")
#
# ax[1][0].set_title("injection long pulse")
# ax[1][0].plot(M1[:,0], M1[:,1], color="k")
#
# ax[0][1].set_title("voltage short pulse")
# ax[0][1].plot(M2[:,0], M2[:,2], color="k")
#
# ax[1][1].set_title("injection short pulse")
# ax[1][1].plot(M2[:,0], M2[:,1], color="k")
#
# ax[0][0].set_ylabel("V (mV)")
# ax[1][0].set_ylabel("I (pA)")
# ax[1][1].set_xlabel("time (ms)")
# ax[1][0].set_xlabel("time (ms)")
# plt.show()

# the parameter to fit Cm, Rm, Ra
CM_IX = 0
RM_IX=1
RA_IX = 2
h.load_file("morphology_correctedZ.hoc")

RM_const = 20000.0
RA_const = 250.0
CM_const = 1.0

print("free params:")

h.attr_praxis(1e-9,1000,0)
opt_vals = h.Vector(3)
opt_vals.x[RM_IX] =22000.0
opt_vals.x[RA_IX] = 200
opt_vals.x[CM_IX] = 1

# allway run the fitting 3 time to avoid stack in local minima
for i in range(3):
    RMSD = h.fit_praxis(efun,opt_vals)

    RM = opt_val=s.x[RM_IX]
    RA = opt_vals.x[RA_IX]
    CM = opt_vals.x[CM_IX]

    print("RMSD", RMSD,", RM",  RM, ", RA",  RA, ", CM",  CM)
    plot_res(CM=CM, RM=RM, RA = RA, save_name="_fit_after_"+str(i+1))

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

