from fit_fun import efun,plot_res,change_model_pas
# from builder import Builder
from open_pickle import read_from_pickle
import numpy as np
from neuron import h,gui
import os
import  matplotlib.pyplot as plt
from simulation import SpineParams
from calculate_F_factor import calculate_F_factor
from add_figure import add_figure
import pickle
import signal
do_calculate_F_factor=True
spine_type="mouse_spine"
initial_folder = "data/fit/"

shrinkage_factor=1.0/0.7
resize_diam_by=1.2

try:os.mkdir(initial_folder)
except FileExistsError:pass
try:os.mkdir(initial_folder+spine_type)
except FileExistsError:pass
# try:os.mkdir(initial_folder+spine_type+"/different_initial_conditions")
# except FileExistsError:pass

if resize_diam_by!=1 and shrinkage_factor!=1:
    initial_folder=initial_folder+spine_type+"/dend*"+str(round(resize_diam_by,2))+' &shrinkage by '+str(round(shrinkage_factor,2))
    try: os.mkdir(initial_folder)
    except FileExistsError: pass
elif resize_diam_by!=1:
    initial_folder=initial_folder+spine_type+"/dend*"+str(round(resize_diam_by,2))
    try: os.mkdir(initial_folder)
    except FileExistsError: pass
elif shrinkage_factor!=1:
    initial_folder=initial_folder+spine_type+"/F_shrinkage="+str(round(shrinkage_factor,2))
    try: os.mkdir(initial_folder)
    except FileExistsError: pass
try:os.mkdir(initial_folder+"/different_initial_conditions")
except FileExistsError:pass
initial_folder=initial_folder+"/different_initial_conditions/"

h.load_file("import3d.hoc")
h.load_file("nrngui.hoc")
h.load_file('stdlib.hoc')
h.load_file("stdgui.hoc")
# h.loadfile("stdrun.hoc")
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

def SIGSEGV_signal_arises(signalNum, stack):
    print(f"{signalNum} : SIGSEGV arises")
    # Your code
signal.signal(signal.SIGSEGV, SIGSEGV_signal_arises)


def change_model_pas(CM=1, RA = 250, RM = 20000.0, E_PAS = -70.0, F_factor = 1.9):
    h.dt = 0.1
    h.distance(0,0.5, sec=soma) # it isn't good beacause it change the synapse distance to the soma
    #h.distance(0, sec=soma)
    for sec in h.allsec(): ##check if we need to insert Ra,cm,g_pas,e_pas to the dendrit or just to the soma
        sec.Ra = RA
        sec.cm = CM*(1.0/0.7)
        sec.g_pas = (1.0 / RM)*(1.0/0.7)
        sec.e_pas = E_PAS
    for sec in cell.dend:
        for seg in sec: #count the number of segment and calclate g_factor and total dend distance,
            # how many segment have diffrent space larger then SPINE_START that decided
            if h.distance(seg) > SPINE_START:
                if do_calculate_F_factor:
                    F_factor=calculate_F_factor(cell,V_head,spine_neck_diam,spine_neck_L)
                seg.cm *= F_factor
                seg.g_pas *= F_factor

def plot_res(RM, RA, CM, save_folder="data/fit/",save_name= "fit",print_full_graph=False):
    change_model_pas(CM=CM, RA=RA, RM=RM, E_PAS = E_PAS,F_factor=F_factor)
    Vvec = h.Vector()
    Tvec = h.Vector()
    Vvec.record(soma(0.5)._ref_v)
    Tvec.record(h._ref_t)
    h.cvode.store_events(Vvec)
    h.run()
    npTvec = np.array(Tvec)
    npVec = np.array(Vvec)
    add_figure("fit "+save_folder.split('/')[-1]+"\nRM="+str(round(RM,1))+",RA="+str(round(RA,1))+",CM="+str(round(CM,2)),'mS','mV')
    plt.plot(npTvec[start_fit:end_fit], npVec[start_fit:end_fit], color = 'r', linestyle ="--")
    plt.plot(T[start_fit:end_fit], V[start_fit:end_fit], color = 'green')
    plt.legend(['NEURON_sim','decay_to_fitting'])
    plt.savefig(save_folder+'/'+save_name+"_decay.png")
    # plt.savefig(save_folder+'/'+save_name+"_decay.pdf")
    plt.close()
    exp_V = V
    npVec = npVec
    npVec = npVec[:len(exp_V)]
    error_1 = np.sqrt(np.sum(np.power(np.mean(exp_V[:2000]) - np.mean(npVec[:2000]), 2)))  # error from mean rest
    error_2 = np.sqrt(np.sum(np.power(exp_V[start_fit+950:end_fit-1000] - npVec[start_fit+950:end_fit-1000], 2)))#/(end_fit-start_fit))  #  error for the decay
    error_3 = np.sqrt(np.sum(np.power(np.mean(exp_V[4100:4900]) - np.mean(npVec[4100:4900]), 2)))  # error for maximal voltage
    error_tot = np.sqrt(np.sum(np.power(exp_V - npVec, 2))/len(exp_V)) # mean square error

    print('error_total=',round(error_tot,3))
    print('error_decay=', round(error_2,3))
    print('error_mean_max_voltage=', round(error_3,3))
    print('error_from_rest=', round(error_1,3))
    return error_2, (error_2 + error_3*10)/960

def efun(vals):
    if RM_IX != -1 :
        if vals.x[RM_IX] > 100000:
            return (1e6)
        RM = vals.x[RM_IX]
    else: RM = RM_const

    if CM_IX != -1:
        if vals.x[CM_IX] >4 :
            return (1e6)
        CM = vals.x[CM_IX]
    else:CM = CM_const

    if RA_IX != -1:
        if vals.x[RA_IX] > 350:
            return (1e6)
        RA = vals.x[RA_IX]
    else:RA = RA_const
    if (CM < 0.3 or RM < 2000 or RA <49):
        return 1e6
    # print('RA:',RA, '   CM:',CM, '   RM:',RM)

    change_model_pas(CM=CM, RA=RA, RM = RM, E_PAS = E_PAS,F_factor=F_factor)
    Vvec = h.Vector()
    Vvec.record(soma(0.5)._ref_v)
    h.run()
    npVec = np.array(Vvec)

    exp_V = V#[int(180.0/h.dt):int(800.0/h.dt)]
    npVec = npVec#[int(180.0/h.dt):int(800.0/h.dt)]
    npVec = npVec[:len(exp_V)]
    error_tot = np.sqrt(np.sum(np.power(exp_V - npVec, 2)))#/len(exp_V)) # mean square error
    error_1 = np.sqrt(np.sum(np.power(np.mean(exp_V[:2000]) - np.mean(npVec[:2000]), 2)))  # error from mean rest
    error_2 = np.sqrt(np.sum(np.power(exp_V[start_fit+950:end_fit-1000] - npVec[start_fit+950:end_fit-1000], 2))) #/(end_fit-start_fit)  #  error for the decay
    error_3 = np.sqrt(np.sum(np.power(np.mean(exp_V[4100:4900]) - np.mean(npVec[4100:4900]), 2)))  # error for maximal voltage

    return error_2 + ((end_fit-1000) - (start_fit+950)) * error_3

def initiate_simulation(self):
    clamp = h.IClamp(self.soma(0.5))  # insert clamp(constant potentientiol) at the soma's center
    clamp.amp = -0.05  ## supopsed to be 0.05nA
    clamp.dur = 200
    clamp.delay = 296
    ######################################################
    # load the data and see what we have
    ######################################################
    V = np.array(self.short_pulse[0])
    T = np.array(self.short_pulse[1])
    T = T - T[0]
    E_PAS = np.mean(V[:2000])
    h.tstop = (T[-1] - T[0]) * 1000
    h.v_init = E_PAS
    h.dt = 0.1
    h.steps_per_ms = h.dt
    return E_PAS, T, V

def fit2short_pulse(cell,short_pulse,folder="",CM=1,RM=10000,RA=100):
    opt_vals = h.Vector(3)
    opt_vals.x[RM_IX] =RM
    opt_vals.x[RA_IX] = RA
    opt_vals.x[CM_IX] = CM
    change_model_pas(CM=CM, RA=RA, RM=RM, E_PAS=E_PAS,F_factor=F_factor)
    plot_res(CM=CM, RM=RM, RA=RA,save_folder=folder, save_name=" before")
    for i in range(3):
        RMSD = h.fit_praxis(efun,opt_vals)   #@# take too much time if the fitting isn't found
        RM = opt_vals.x[RM_IX]
        RA = opt_vals.x[RA_IX]
        CM = opt_vals.x[CM_IX]

        print("RMSD", RMSD,", RM",  RM, ", RA",  RA, ", CM",  CM)
        error2,precent_error=plot_res(CM=CM, RM=RM, RA=RA, save_folder=folder,save_name="_fit_after_" + str(i + 1))
    pickle.dump({
        "RM": RM,
        "RA": RA,
        "CM": CM,
        "error":[RMSD,error2,precent_error]
    }, open(folder+"/fit_result.p", "wb"))
    return {"CM": CM,"RM": RM,"RA": RA,"error":[error2,precent_error,RMSD]}
if __name__=='__main__':
    I = -50
    cell_file= "05_08_A_01062017_Splice_shrink_FINISHED_LABEL_Bluecell_spinec91.ASC"
    short_pulse_file="data/short_pulse/mean_short_pulse_with_parameters.p"
    short_pulse=read_from_pickle(short_pulse_file)
    class Cell:pass
    cell = mkcell(cell_file)
    sp = h.PlotShape()
    sp.show(0)  # show diameters
    for sec in h.allsec():
        sec.diam = sec.diam*resize_diam_by
        sec.L*=shrinkage_factor
    ## delete all the axons
    for sec in cell.axon:
        h.delete_section(sec=sec)
    for sec in h.allsec():
        sec.insert('pas') # insert passive property
        sec.nseg = int(sec.L/10)+1  #decide that the number of segment will be 21 with the same distances

    short_pulse=read_from_pickle(short_pulse_file)

    SPINE_START = 60
    start_fit =  2000
    end_fit =  4900
    spine=SpineParams()
    F_factor = 2.03#calculate_F_factor(cell, spine.V_head, spine.neck_diam, spine.neck_L)
    soma=cell.soma[0]
    clamp = h.IClamp(soma(0.5)) # insert clamp(constant potentientiol) at the soma's center
    clamp.amp = -0.05
    clamp.dur = 200
    clamp.delay = 296
    ######################################################
    # load the data and see what we have
    ######################################################
    V = np.array(short_pulse['mean'][0])
    T = np.array(short_pulse['mean'][1].rescale('ms'))
    T = T-T[0]

    E_PAS = short_pulse['E_pas']
    h.tstop = (T[-1]-T[0])
    h.v_init=E_PAS
    h.dt = 0.1
    h.steps_per_ms = h.dt
    CM_IX = 2
    RM_IX=0
    RA_IX = 1

    RM_const = 60000.0
    RA_const = 100.0
    CM_const = 1.0
    print("free params:")
    h.attr_praxis(1e-9,1000,0)

    CM = 1  # 2/2
    RM = 5684*2  # *2
    RA = 100

    ra_folder = initial_folder + "/RA0_50:100:0.5"
    try:os.mkdir(ra_folder)
    except FileExistsError:pass
    RAs = np.arange(50,100,0.5)
    solution_RA0={}
    for ra in RAs:
        folder = ra_folder + "/RA0=" + str(ra)
        try:os.mkdir(folder)
        except FileExistsError:pass
        solution_RA0["RA0=" + str(ra)] = fit2short_pulse(cell, short_pulse, folder=folder, CM=CM, RM=RM, RA=ra)
        pickle.dump(solution_RA0, open(ra_folder + "/RA0_fit_results.p", "wb"))

    ra_folder = initial_folder + "/RA0_100:300:2"
    try:os.mkdir(ra_folder)
    except FileExistsError:pass
    RAs = np.arange(100,300,2.)
    solution_RA0={}
    for ra in RAs:
        folder = ra_folder + "/RA0=" + str(ra)
        try:os.mkdir(folder)
        except FileExistsError:pass
        solution_RA0["RA0=" + str(ra)] = fit2short_pulse(cell, short_pulse, folder=folder, CM=CM, RM=RM, RA=ra)
    pickle.dump(solution_RA0, open(ra_folder + "/RA0_fit_results.p", "wb"))
    # cm_folder = initial_folder+"/CM0"
    # try:os.mkdir(cm_folder)
    # except FileExistsError:pass
    # CMs=[0.5,0.8,1,1.2,1.4,1.8,2,2.5,3]
    # solution_CM0={}
    # for cm in tqdm(CMs):
    #     folder = cm_folder+"/CM0="+str(cm)
    #     try:os.mkdir(folder)
    #     except FileExistsError:pass
    #     solution_CM0["CM0="+str(cm)]=fit2short_pulse(cell,short_pulse,folder=folder,CM=cm,RM=RM,RA=RA)
    # pickle.dump(solution_CM0 , open(cm_folder+"/CM0_fit_results.p", "wb"))

    # rm_folder = initial_folder+"/RM0"
    # try:os.mkdir(rm_folder)
    # except FileExistsError:pass
    # RMs=[5000,10000,15000,20000,25000,30000,50000,80000]
    # solution_RM0={}
    # for rm in RMs:
    #     folder = rm_folder+"/RM0="+str(rm)
    #     try:os.mkdir(folder)
    #     except FileExistsError:pass
    #     solution_RM0["RM0="+str(rm)]=fit2short_pulse(cell,short_pulse,folder=folder,CM=CM,RM=rm,RA=RA)
    # pickle.dump(solution_RM0 , open(rm_folder+"/RM0_fit_results.p", "wb"))





    # CM_RA_RM_folder = initial_folder + "/CM0_RM0_RA0"
    # try:os.mkdir(CM_RA_RM_folder)
    # except FileExistsError:pass
    # solution={}
    # for cm in CMs:
    #     for rm in RMs:
    #         for ra in RAs:
    #             folder = CM_RA_RM_folder + "/CM0="+str(cm)+" RM0="+str(rm)+" RA0=" + str(ra)
    #             try:os.mkdir(folder)
    #             except FileExistsError:pass
    #             solution["CM0="+str(cm)+" RM0="+str(rm)+" RA0=" + str(ra)] = fit2short_pulse(cell, short_pulse, folder=folder, CM=cm, RM=rm, RA=ra)
    # pickle.dump(solution, open(CM_RA_RM_folder + "/CM_RM_RA_fit_results.p", "wb"))
