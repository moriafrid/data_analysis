from open_MOO_after_fit import OPEN_RES
import numpy as np
# from neuron import h
import matplotlib.pyplot as plt
import os
from glob import glob
folder_='AMPA&NMDA/'
try:os.mkdir(folder_)
except FileExistsError:pass
for type in ['hand','outomatic']:
    base='../test_moo/groger_spine/MOO_opt_folder_same_w_peel_syn_spines_featues_cm_g_pas_first_seed_123456/05_08_A_01062017_fit.py_solution_'+type
    if type=='outomatic':
        base=base+'/run_78_'
    elif type=='hand':
        base=base+'/run_75_'
    for model_place in glob(base+'*'):
        loader = OPEN_RES(res_pos=model_place+'/')

        model=loader.get_model()
        h=loader.sim.neuron.h
        netstim = h.NetStim()  # the location of the NetStim does not matter
        netstim.number = 1
        netstim.start = 200
        netstim.noise = 0
        spine, syn_obj = loader.create_synapse(model.dend[82], 0.165, netstim=netstim)
        h.tstop = 400
        time = h.Vector()
        time.record(h._ref_t)
        V_soma = h.Vector()
        V_soma.record(model.soma[0](0.5)._ref_v)
        h.dt = 0.1
        h.steps_per_ms = 1.0/h.dt
        h.run()
        V_soma_All = np.array(V_soma)[1900:]
        time_all = np.array(time)[1900:]

        syn_obj[1][1].weight[0]=0
        h.dt=0.1
        h.steps_per_ms = 1.0/h.dt
        h.run()
        V_soma_AMPA = np.array(V_soma)[1900:]
        time_AMPA = np.array(time)[1900:]
        V_NMDA = V_soma_All-V_soma_AMPA
        from add_figure import add_figure
        passive_propert_title='Rm='+str(round(1.0/model.dend[82].g_pas,2)) +' Ra='+str(round(model.dend[82].Ra,2))+' Cm='+str(round(model.dend[82].cm,2))
        add_figure('AMPA and NMDA impact on voltage '+" ".join(model_place.split('/')[-1].split('_')[2:])+'\n'+passive_propert_title,'time[ms]','Voltage[mV]')
        plt.plot(time_all, V_soma_All, color='g', lw=5,label='all',alpha=0.4)
        plt.plot(time_all, V_soma_AMPA, color='b', lw=2,linestyle='--', label='AMPA',alpha=0.8)
        plt.plot(time_all, V_NMDA+V_soma_All[0],lw=2, color='r', linestyle='--', label='NMDA',alpha=0.8)
        plt.legend()
        plt.savefig(folder_+model_place.split('/')[-1]+'.png')
