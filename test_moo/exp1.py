from open_MOO_after_fit import OPEN_RES
import numpy as np
from neuron import h
import matplotlib.pyplot as plt
import os
from glob import glob
folder_='Voltage Spine&Soma/'
try:os.mkdir(folder_)
except FileExistsError:pass
for type in ['hand','outomatic']:
    base='../test_moo/groger_spine/MOO_opt_folder_same_w_peel_syn_spines_featues_cm_g_pas_first_seed_123456/05_08_A_01062017_fit.py_solution_'+type
    if type=='outomatic':
        base=base+'/run_78_'
    elif type=='hand':
        base=base+'/run_76_'
    for model_place in glob(base+'*'):
        loader = OPEN_RES(res_pos=model_place+'/')
        model=loader.get_model()
        netstim = h.NetStim()  # the location of the NetStim does not matter
        netstim.number = 1
        netstim.start = 200
        netstim.noise = 0
        spine, syn_obj = loader.create_synapse(model.dend[82], 0.165, netstim=netstim)
        h.tstop = 300
        time = h.Vector()
        time.record(h._ref_t)
        V_soma = h.Vector()
        V_soma.record(model.soma[0](0.5)._ref_v)
        V_spine = h.Vector()
        V_spine.record(spine[1](1)._ref_v)
        h.run()
        from add_figure import add_figure
        passive_propert_title='Rm='+str(round(1.0/model.dend[82].g_pas,2)) +' Ra='+str(round(model.dend[82].Ra,2))+' Cm='+str(round(model.dend[82].cm,2))
        V_soma=np.array(V_soma)[1900:]
        V_spine=np.array(V_spine)[1900:]
        time=np.array(time)[1900:]
        add_figure('Voltage in Spine and Soma '+" ".join(model_place.split('/')[-1].split('_')[2:])+'\n'+passive_propert_title,'mV','mS')
        plt.plot(time, V_soma,label='V_soma max:'+str(round(np.amax(V_soma[1000:]),2)))
        plt.plot(time, V_spine,label='V_spine max:'+str(round(np.amax(V_spine[1000:]),2)))
        plt.legend()
        plt.savefig(folder_+model_place.split('/')[-1]+'.png')
