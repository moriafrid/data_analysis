from open_MOO_after_fit import OPEN_RES
import numpy as np
from neuron import h
import matplotlib.pyplot as plt
import os
from glob import glob
from tqdm import tqdm
from add_figure import add_figure

folder_='hall_of_fame_together/'
type='outomatic'
try:os.mkdir(folder_)
except FileExistsError:pass
for type in ['outomatic']:
    base='../test_moo/groger_spine/MOO_opt_folder_same_w_peel_syn_spines_featues_cm_g_pas_first_seed_123456/05_08_A_01062017_fit.py_solution_'+type
    if type=='outomatic':
        base=base+'/run_78_'
    elif type=='hand':
        base=base+'/run_75_'
    for model_place in glob(base+'*'):
        add_figure('AMPA and NMDA impact on voltage ' + model_place.split('/')[-1].split('_')[-1],'mV', 'mS')
        print(model_place.split('/')[-2])
        for i in tqdm(range(10)):
            loader = OPEN_RES(res_pos=model_place + '/')
            model=loader.get_model()
            netstim = h.NetStim()  # the location of the NetStim does not matter
            netstim.number = 1
            netstim.start = 200
            netstim.noise = 0
            spine, syn_obj = loader.create_synapse(model.dend[82], 0.165, netstim=netstim, hall_of_fame_num=i)
            h.tstop = 300
            time = h.Vector()
            time.record(h._ref_t)
            V_spine = h.Vector()
            V_spine.record(spine[1](1)._ref_v)
            h.dt = 0.1
            h.steps_per_ms = 1.0 / h.dt
            h.run()
            V_spine=np.array(V_spine)[1900:]
            time=np.array(time)[1900:]
            plt.plot(time, V_spine, label='hall_'+str(i),alpha=0.1)
        passive_propert_title='Rm='+str(round(1.0/model.dend[82].g_pas,2)) +' Ra='+str(round(model.dend[82].Ra,2))+' Cm='+str(round(model.dend[82].cm,2))
        plt.title('AMPA and NMDA impact on voltage ' +" ".join(model_place.split('/')[-1].split('_')[2:]) + '\n' + passive_propert_title)
        plt.legend()
        plt.savefig(folder_+model_place.split('/')[-1]+'.png')
