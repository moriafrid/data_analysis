import numpy as np
from neuron import h,gui
from builder import Builder
import matplotlib.pyplot as plt
from open_pickle import read_from_pickle
class FitModel:
    def __init__(self,
                 cell_file,
                 short_pulse_file,
                 V_head = 0.14,
                 spine_neck_diam = 0.164,
                 spine_neck_L = 0.782,
                 SPINE_START = 60,
                 start_fit = 0,
                 end_fit = 4900):
        self.cell = Builder.initiete_cell(cell_file)
        self.short_pulse = read_from_pickle(short_pulse_file)
        self.soma = self.cell.soma[0]


    def change_model_pas(self,CM=1, RA = 250, RM = 20000.0, E_PAS = -70.0, F_factor = {}):
        h.dt = 0.1
        h.distance(0,0.5, sec=self.soma) # it isn't good beacause it change the synapse distance to the soma
        #h.distance(0, sec=soma)
        for sec in h.allsec(): ##check if we need to insert Ra,cm,g_pas,e_pas to the dendrit or just to the soma
       #@# cell.axon is not exist in hoc object
            sec.Ra = RA
            sec.cm = CM
            sec.g_pas = 1.0 / RM
            sec.e_pas = E_PAS
        for sec in self.cell.dend:
            for seg in sec: #count the number of segment and calclate g_factor and total dend distance,
                # how many segment have diffrent space larger then SPINE_START that decided
                if h.distance(seg) > SPINE_START:
                    if do_calculate_F_factor:
                        F_factor=calculate_F_factor(self.cell,V_head,spine_neck_diam,spine_neck_L)
                    else:
                        F_factor = 1.9#2.0 # F_factor[sec]
                    seg.cm *= F_factor
                    seg.g_pas *= F_factor
    ## e_pas is the equilibrium potential of the passive current

    def plot_res(self,RM, RA, CM, save_folder="data/fit/",save_name= "fit",print_full_graph=False):
        FitModel.change_model_pas(CM=CM, RA=RA, RM=RM, E_PAS = self.E_PAS)
        Vvec = h.Vector() #cerat vector to record on it
        Tvec = h.Vector() #cerat vector to record on it
        Vvec.record(self.soma(0.5)._ref_v) #where to recprd
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
        return save_folder

    def efun(self,vals):
        if self.RM_IX != -1 :
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
        if (CM < 0.3 or RM < 2000 or RA <1):
            return 1e6
        # print('RA:',RA, '   CM:',CM, '   RM:',RM)

        FitModel.change_model_pas(CM=CM, RA=RA, RM = RM, E_PAS = self.E_PAS)
        Vvec = h.Vector()
        Vvec.record(self.soma(0.5)._ref_v)
        h.run()
        npVec = np.array(Vvec)

        exp_V = V#[int(180.0/h.dt):int(800.0/h.dt)]
        npVec = npVec#[int(180.0/h.dt):int(800.0/h.dt)]
        npVec = npVec[:len(exp_V)]
        error_tot = np.sqrt(np.sum(np.power(exp_V - npVec, 2)))#/len(exp_V)) # mean square error
        error_1 = np.sqrt(np.sum(np.power(np.mean(exp_V[:2000]) - np.mean(npVec[:2000]), 2)))  # error from mean rest
        error_2 = np.sqrt(np.sum(np.power(exp_V[start_fit:end_fit] - npVec[start_fit:end_fit], 2))) #/(end_fit-start_fit)  #  error for the decay
        error_3 = np.sqrt(np.sum(np.power(np.mean(exp_V[4100:4900]) - np.mean(npVec[4100:4900]), 2)))  # error for maximal voltage

        return error_2 + (end_fit-start_fit)*error_3

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


class PassiveParams:
    RM_const = 60000.0
    RA_const = 250.0
    CM_const = 1.0

class SpineParams:
    V_head = 0.14
    neck_diam = 0.164
    neck_L = 0.782


class Simulation_parameters:
    SPINE_START = 60
    start_fit = 0  # 2960
    end_fit = 4900  # 3960
