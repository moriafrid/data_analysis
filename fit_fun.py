import signal
from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt
from open_pickle import read_from_pickle
import os
from add_figure import add_figure
import pickle

def change_model_pas(CM=1, RA = 250, RM = 20000.0, E_PAS = -70.0, F_factor = {}):
    #input the neuron property    h.dt = 0.1

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
                else:
                    F_factor = 1.9#2.0 # F_factor[sec]
                seg.cm *= F_factor
                seg.g_pas *= F_factor
def plot_res(RM, RA, CM, save_name= "fit",print_full_graph=False):
    folder_="data/fit/"
    try: os.mkdir(folder_)
    except FileExistsError: pass
    if resize_diam_by==1:
        try: os.mkdir(folder_+str(I)+'pA/')
        except FileExistsError: pass
        save_folder=folder_+str(I)+'pA/'
    else:
        try: os.mkdir(folder_+'dend*'+str(resize_diam_by)+'_'+str(I)+'pA/')
        except FileExistsError: pass
        save_folder=folder_+'dend*'+str(resize_diam_by)+'_'+str(I)+'pA/'
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
        if vals.x[CM_IX] >4 :
            return (1e6)
        CM = vals.x[CM_IX]
    else:CM = CM_const

    if RA_IX != -1:
        if vals.x[RA_IX] > 300:
            return (1e6)
        RA = vals.x[RA_IX]
    else:RA = RA_const
    if (CM < 0.3 or RM < 2000 or RA <1):
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
def fit2short_pulse(cell,short_pulse,folder="",CM=1,RM=10000,RA=100):
    opt_vals = h.Vector(3)
    opt_vals.x[RM_IX] =RM
    opt_vals.x[RA_IX] = RA
    opt_vals.x[CM_IX] = CM
    change_model_pas(CM=CM, RA=RA, RM=RM, E_PAS=E_PAS)

    imp = h.Impedance(sec=soma)
    imp.loc(soma(0.5))
    imp.compute(0)
    imp.input(0)

    plot_res(CM=CM, RM=RM, RA=RA, save_name="before")
    print('the initial impadence is', imp.input(0))
    # allway run the fitting 3 time to avoid stack in local minima
    for i in range(3):
        RMSD = h.fit_praxis(efun,opt_vals)   #@# take too much time if the fitting isn't found
        RM = opt_vals.x[RM_IX]
        RA = opt_vals.x[RA_IX]
        CM = opt_vals.x[CM_IX]

        print("RMSD", RMSD,", RM",  RM, ", RA",  RA, ", CM",  CM)
        plot_res(CM=CM, RM=RM, RA=RA, save_name="_fit_after_" + str(i + 1))
        imp.compute(0)
        print('the impadence is',imp.input(0))
        pickle.dump({
            "RM": RM,
            "RA": RA,
            "CM": CM
        }, open(folder + '/' + "final_result_dend*" + str(resize_diam_by) + ".p", "wb"))
        return {"RMSD": RMSD, "RM":RM, "RA":RA, "CM":  CM}

if __name__=='__main__':
    folder_="data/fit/"
    try: os.mkdir(folder_)
    except FileExistsError: pass
    SPINE_START = 60
    resize_diam_by = 1
    hz = 10000
    start_fit = 0  # 2960
    end_fit = 4900  # 3960
    from simulation import SpineParams
    from builder import Builder
    from calculate_F_factor import calculate_F_factor
    cell_file = "05_08_A_01062017_Splice_shrink_FINISHED_LABEL_Bluecell_spinec91.ASC"
    cell=Builder.initiete_cell(cell_file)
    short_pulse_file="data/short_pulse/clear_short_pulse.p"
    short_pulse=read_from_pickle(short_pulse_file)
    soma=cell.soma[0]
    spine=SpineParams()
    F_factor = calculate_F_factor(cell, spine.V_head, spine.spine_neck_diam, spine.spine_neck_L)

    clamp = h.IClamp(soma(0.5)) # insert clamp(constant potentientiol) at the soma's center
    clamp.amp = -0.05 ## supopsed to be 0.05nA
    clamp.dur = 200
    clamp.delay = 296
    ######################################################
    # load the data and see what we have
    ######################################################
    V = np.array(short_pulse[0])
    T = np.array(short_pulse[21])
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

    solution = fit2short_pulse(cell, short_pulse, folder=folder, RA=RA,CM=CM,RM=RM)


