from fit_fun import efun,plot_res,change_model_pas
from builder import Builder
from open_pickle import read_from_pickle
import numpy as np
from neuron import h,gui
import os
import  matplotlib.pyplot as plt
from simulation import SpineParams
from calculate_F_factor import calculate_F_factor
from add_figure import add_figure
import pickle
from math import pi
initial_folder = "data/fit/"
do_calculate_F_factor=True
spine_type="mouse_spine"

shrinkage_factor=1.2#1.0/0.7
resize_diam_by=1


try:os.mkdir(initial_folder+spine_type)
except FileExistsError:pass

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
else:
    initial_folder=initial_folder+spine_type+"/no change"
    try: os.mkdir(initial_folder)
    except FileExistsError: pass

try:os.mkdir(initial_folder+"/const_param")
except FileExistsError:pass
try:os.mkdir(initial_folder+"/const_param/RA")
except FileExistsError:pass
initial_folder = initial_folder+"/const_param/RA"

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


def change_model_pas(CM=1, RA = 250, RM = 20000.0, E_PAS = -70.0, F_factor = 1.9):
    h.dt = 0.1
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
                seg.cm *= F_factor
                seg.g_pas *= F_factor

def plot_res(RM, RA, CM, save_folder="data/fit/",save_name= "fit"):
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
    plt.plot(npTvec[start_fit:end_fit], npVec[start_fit:end_fit], color = 'r', linestyle ="--",alpha=0.3)
    plt.plot(npTvec[start_fit+950:end_fit-1000], npVec[start_fit+950:end_fit-1000], color = 'b',alpha=0.3)
    plt.plot(T[start_fit:end_fit], V[start_fit:end_fit], color = 'green',alpha=0.3)
    plt.legend(['NEURON_sim','decay_to_fitting'])
    plt.savefig(save_folder+'/'+save_name+"_decay.png")
    # plt.savefig(save_folder+'/'+save_name+"_decay.pdf")
    plt.close()
    exp_V = V
    npVec = npVec
    npVec = npVec[:len(exp_V)]
    error_1 = np.sqrt(np.sum(np.power(np.mean(exp_V[:2000]) - np.mean(npVec[:2000]), 2)))  # error from mean rest
    error_2 = np.sqrt(np.sum(np.power(exp_V[start_fit+950:end_fit-1000] - npVec[start_fit+950:end_fit-1000], 2)))#/((end_fit-1000) - (start_fit+950)))  #  error for the decay
    error_3 = np.sqrt(np.sum(np.power(np.mean(exp_V[4100:4900]) - np.mean(npVec[4100:4900]), 2)))  # error for maximal voltage
    error_tot = np.sqrt(np.sum(np.power(exp_V - npVec, 2))/len(exp_V)) # mean square error

    print('error_total=',round(error_tot,3))
    print('error_decay=', round(error_2,3))
    print('error_mean_max_voltage=', round(error_3,3))
    print('error_from_rest=', round(error_1,3))
    return error_2, (error_2 + error_3*10)/960
def errors_Rinput(RM,RA,CM,E_PAS):
    change_model_pas(CM=CM, RA=RA, RM=RM, E_PAS = E_PAS,F_factor=F_factor)
    Vvec = h.Vector()
    Tvec = h.Vector()
    Vvec.record(soma(0.5)._ref_v)
    Tvec.record(h._ref_t)
    h.cvode.store_events(Vvec)
    h.dt=0.1
    h.run()
    npTvec = np.array(Tvec)
    npVec = np.array(Vvec)
    exp_V = V
    npVec = npVec
    npVec = npVec[:len(exp_V)]
    error_3 = np.sqrt(np.sum(np.power(np.mean(exp_V[4100:4900]) - np.mean(npVec[4100:4900]), 2)))  # error for maximal voltage
    # print('error_mean_max_voltage=', round(error_3,3))
    return error_3

if __name__=='__main__':
    I = -50
    cell_file= "05_08_A_01062017_Splice_shrink_FINISHED_LABEL_Bluecell_spinec91.ASC"
    short_pulse_file="data/short_pulse/mean_short_pulse_with_parameters.p"
    short_pulse=read_from_pickle(short_pulse_file)
    cell=Builder.initiate_cell(cell_file)
    for sec in h.allsec():
        sec.diam = sec.diam*resize_diam_by
        sec.L*=shrinkage_factor
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
    imp = h.Impedance(sec=soma)
    imp.loc(soma(0.5))
    RA=np.arange(200, 300, 2)
    tau_m=25716#ms*10^3=micro s
    d=soma.diam
    ra_error=[]
    params_dict=[]
    precent_erors=[]
    ra_error_next=[]
    for ra in RA:
        RM=6000
        CM=tau_m/RM
        imp.compute(0)
        Rin = imp.input(0)
        error_last=errors_Rinput(RM, ra, CM,E_PAS)
        error_next=error_last
        while error_next<=error_last:
            RM+=20
            CM=tau_m/RM
            change_model_pas(CM=CM, RA = ra, RM = RM, E_PAS = E_PAS, F_factor = F_factor)
            imp.compute(0)
            Rin=imp.input(0)
            # print('Rin='+str(round(rin,3)))
            # RM=(Rin*pi)**2/4*d**3*ra
            error_last=error_next
            error_next=errors_Rinput(RM, ra, CM,E_PAS)
        print(RM)
        RM-=35
        CM = tau_m / RM
        change_model_pas(CM=CM, RA=ra, RM=RM, E_PAS=E_PAS, F_factor=F_factor)
        imp.compute(0)
        Rin = imp.input(0)
        error_last = errors_Rinput(RM, ra, CM, E_PAS)
        error_next = errors_Rinput(RM, ra, CM, E_PAS)
        while error_next<=error_last:
            RM+=1
            CM=tau_m/RM
            change_model_pas(CM=CM, RA = ra, RM = RM, E_PAS = E_PAS, F_factor = F_factor)
            imp.compute(0)
            Rin=imp.input(0)
            # print('Rin='+str(round(rin,3)))
            # RM=(Rin*pi)**2/4*d**3*ra
            error_last=error_next
            error_next=errors_Rinput(RM, ra, CM,E_PAS)
        print('Rinput for Ra='+str(ra)+' is '+str(round(Rin,2)))
        ra_error2,precent_eror=plot_res(RM, ra, CM, save_folder=initial_folder, save_name="fit for RA=" + str(round(ra, 2)))
        ra_error.append(ra_error2)
        precent_erors.append(precent_eror)
        ra_error_next.append(error_next)
        params_dict.append({'RM': RM, 'RA': ra, 'CM': CM})
    pickle.dump({'RA':RA,'errors':[ra_error,precent_erors],'params':params_dict}, open(initial_folder + '/Ra_const_errors200:300.p', "wb"))
    add_figure('RA_errors','RA','errors')
    plt.plot(RA,ra_error)
    plt.savefig(initial_folder + '/Ra_const_errors1.png')
    plt.add_figure('RA_errors','RA','ra_next_eror')
    plt.plot(RA,ra_error_next)
    plt.savefig(initial_folder + '/Ra_const_errors2.png')

