from open_MOO_after_fit import OPEN_RES
import numpy as np
# from neuron import h
import matplotlib.pyplot as plt
import os
from glob import glob
folder_='g_max/'
type='outomatic'#'hand'
try:os.mkdir(folder_)
except FileExistsError:pass
for type in ['hand','outomatic']:
    base='../test_moo/groger_spine/MOO_opt_folder_same_w_peel_syn_spines_featues_cm_g_pas_first_seed_123456/05_08_A_01062017_fit.py_solution_'+type
    if type=='outomatic':
        base=base+'/run_78_'
    elif type=='hand':
        base=base+'/run_75_'
    for model_place in glob(base + '*'):
        loader = OPEN_RES(res_pos=model_place + '/')

        model=loader.get_model()
        h=loader.sim.neuron.h
        netstim = h.NetStim()  # the location of the NetStim does not matter
        netstim.number = 1
        netstim.start = 200
        netstim.noise = 0
        spine, syn_obj = loader.create_synapse(model.dend[82], 0.165, netstim=netstim)
        h.tstop = 500
        time = h.Vector()
        time.record(h._ref_t)
        V_spine = h.Vector()
        V_spine.record(spine[1](1)._ref_v)
        g_spine_NMDA = h.Vector()
        g_spine_NMDA.record(syn_obj[1][0]._ref_g_NMDA)
        # g_spine_AMPA = h.Vector()
        # g_spine_AMPA.record(syn_obj[1][0]._ref_g_AMPA)
        h.dt = 0.1
        h.steps_per_ms = 1.0/h.dt
        h.run()
        V_spine_All = np.array(V_spine)[1900:]
        g_spine_All = np.array(g_spine_NMDA)[1900:]
        # g_spine_All = np.array(g_spine_AMPA)[1900:]
        time_all = np.array(time)[1900:]
        from add_figure import add_figure

        passive_propert_title = 'Rm=' + str(round(1.0 / model.dend[82].g_pas, 2)) + ' Ra=' + str(round(model.dend[82].Ra, 2)) + ' Cm=' + str(round(model.dend[82].cm, 2))
        add_figure('NMDA g ' +" ".join(model_place.split('/')[-1].split('_')[2:]) + '\n' + passive_propert_title,'time [s]', 'leakness [Simans]')
        # plt.plot(time_all, V_spine_All-V_spine_All[0], color='k', label='all',alpha=0.3)
        plt.plot(time_all, g_spine_All, color='red', linestyle='--', label='NMDA gmax')
        # plt.plot(time_all, g_spine_All, color='b', linestyle='--', label='AMPA gmax')
        plt.legend()
        plt.savefig(folder_+model_place.split('/')[-1]+'.png')
