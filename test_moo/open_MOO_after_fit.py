#! /ems/elsc-labs/segev-i/yoni.leibner/anaconda2/bin/ipython
#70
from __future__ import print_function
import bluepyopt as bpopt
import bluepyopt.ephys as ephys
import pprint
import numpy as np
import os, pickle
from add_figure import add_figure
from open_pickle import read_from_pickle
from calculate_F_factor import calculate_F_factor
import signal
def SIGSEGV_signal_arises(signalNum, stack):
    print(f"{signalNum} : SIGSEGV arises")
signal.signal(signal.SIGSEGV, SIGSEGV_signal_arises)


class OPEN_RES():
    def __init__(self, res_pos='/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/test_moo/groger_spine/MOO_opt_folder_same_w_peel_syn_spines_featues_cm_g_pas_first_seed_123456/05_08_A_01062017/run_74/'):
        from extraClasses import NrnSegmentSomaDistanceScaler_, NrnSectionParameterPas, neuron_start_time, \
            EFeatureImpadance, EFeaturePeak, EFeaturePeakTime, EFeatureRDSM, NrnNetstimWeightParameter, \
            SweepProtocolRin2
        self.res_pos = res_pos
        self.hall = hall = pickle.load(open(self.res_pos + 'hall_of_fame.p', 'rb'))
        self.morph_path = hall['model'].split('morphology')[1].split('\n')[1].strip() # remove spaces
        self.fixed_params_res = dict()
        self.optimization_params_res = dict()
        for line in hall['model'].split('params:')[1].split('\n'):
            param_name = line.split(':')[0]
            if len(line.split('='))>=2 and line.split('=')[1].find('[')==-1:
                self.fixed_params_res[param_name.strip()] = float(line.split('=')[1])
        hall_of_phase_results = np.array(hall['hall_of_fame'].items)
        for i, param in enumerate(hall['parameters']):
            self.optimization_params_res[param] = hall_of_phase_results[:, i]

        self.mechanisms={}
        somatic_loc = ephys.locations.NrnSeclistLocation('somatic', seclist_name='somatic')
        basal_loc = ephys.locations.NrnSeclistLocation('basal', seclist_name='basal')
        apical_loc = ephys.locations.NrnSeclistLocation('apical', seclist_name='apical')
        axonal_loc = ephys.locations.NrnSeclistLocation('axonal', seclist_name='axonal')
        location_dict = {'all': [somatic_loc, basal_loc, apical_loc, axonal_loc],
                         'somatic': [somatic_loc],
                         'apical': [apical_loc],
                         'basal': [basal_loc],
                         'axonal': [axonal_loc]}

        # self.spine_properties=hall['spine']
        self.morphology = ephys.morphologies.NrnFileMorphology(self.morph_path, do_replace_axon=True,
                                                      do_resize_dend=False,nseg_frequency=40)
        self.mechanism_list=[]
        sec_list=location_dict ['all']
        # for mech in self.mechanisms.keys():
        self.mechanism_list.append(ephys.mechanisms.NrnMODMechanism(name='pas', prefix='pas', locations=sec_list))
        F_FACTOR_DISTANCE = NrnSegmentSomaDistanceScaler_(name='spine_factor',
                                                          dist_thresh_apical=60,
                                                          dist_thresh_basal=60,
                                                          F_factor=hall['spine']['F_factor'])
        self.parameters_list=[]
        for parameter in self.fixed_params_res.keys():
            if parameter in ['cm', 'g_pas']:
                self.parameters_list.append(ephys.parameters.NrnSectionParameter(name=parameter, param_name=parameter,value_scaler=F_FACTOR_DISTANCE, value=self.fixed_params_res[parameter], locations=sec_list,frozen=True))
            elif parameter in ['Ra', 'e_pas']:
                self.parameters_list.append(ephys.parameters.NrnSectionParameter(name=parameter, param_name=parameter, value=self.fixed_params_res[parameter], locations=sec_list,frozen=True))
        self.sim=ephys.simulators.NrnSimulator(cvode_active=False)
        self.model = ephys.models.CellModel('Model', morph=self.morphology, mechs=self.mechanism_list,
                                   params=self.parameters_list)
        self.model.instantiate(self.sim)
        self.hoc_model = self.sim.neuron.h.Model[-1]

    def get_model(self):
        return self.hoc_model

    def create_synapse(self, sec, pos, number = 0,
                       neck_diam = 0.25, neck_length = 1.35,
                       head_diam = 0.944, hall_of_fame_num = 0, netstim=None):

        spine = self.create_spine(sec, pos, number = number,
                       neck_diam = neck_diam, neck_length = neck_length,
                       head_diam = head_diam)
        syn_obj = self._add_syn_on_sec(spine[1], 1, hall_of_fame_num=hall_of_fame_num, netstim=netstim)
        return spine, syn_obj


    def create_spine(self, sec, pos, number=0, neck_diam=0.25, neck_length=1.35,
                     head_diam=0.944):  # np.sqrt(2.8/(4*np.pi))
        neck = self.sim.neuron.h.Section(name="spineNeck" + str(number))
        head = self.sim.neuron.h.Section(name="spineHead" + str(number))
        self.hoc_model.all.append(neck)
        self.hoc_model.all.append(head)
        neck.L = neck_length
        neck.diam = neck_diam
        head.L = head.diam = head_diam
        head.connect(neck(1))
        neck.connect(sec(pos))
        self.sim.neuron.h("access " + str(neck.hoc_internal_name()))
        self.hoc_model.all.append(neck)
        self.sim.neuron.h.pop_section()
        self.sim.neuron.h("access " + str(head.hoc_internal_name()))
        self.hoc_model.all.append(head)
        self.sim.neuron.h.pop_section()
        for sec in [neck, head]:
            sec.insert("pas")
            sec.g_pas = self.get_param('g_pas')
            sec.cm = self.get_param('cm')
            sec.Ra = self.get_param('Ra')
        return [neck, head]

    def get_param(self, param_name, hall_num=0):
        if param_name in self.fixed_params_res.keys():
            return self.fixed_params_res[param_name]
        elif param_name in self.optimization_params_res.keys():
            return self.optimization_params_res[param_name][hall_num]
        raise Exception('parameter name isnt correct: '+str(param_name))

    def _add_syn_on_sec(self, sec, pos=1, netstim=None, hall_of_fame_num=0):
        if netstim == None:
            raise Exception('we need netstim!:)')
        AMPA_PART = self.sim.neuron.h.Exp2Syn(sec(pos))
        NMDA_PART = self.sim.neuron.h.NMDA(sec(pos))
        AMPA_PART.tau1 = self.get_param('exp2syn_tau1', hall_of_fame_num)
        AMPA_PART.tau2 = self.get_param('exp2syn_tau2', hall_of_fame_num)
        AMPA_PART.e = 0
        NMDA_PART.e = 0
        NMDA_PART.tau_r_NMDA=self.get_param('NMDA_tau_r_NMDA', hall_of_fame_num)
        NMDA_PART.tau_d_NMDA=self.get_param('NMDA_tau_d_NMDA', hall_of_fame_num)
        NMDA_PART.n_NMDA=self.get_param('NMDA_n_NMDA', hall_of_fame_num)
        NMDA_PART.gama_NMDA=self.get_param('NMDA_gama_NMDA', hall_of_fame_num)
        netcon_AMPA = self.sim.neuron.h.NetCon(netstim, AMPA_PART)
        netcon_NMDA = self.sim.neuron.h.NetCon(netstim, NMDA_PART)
        netcon_AMPA.weight[0] = self.get_param('weight_AMPA', hall_of_fame_num)
        netcon_NMDA.weight[0] = self.get_param('weight_NMDA', hall_of_fame_num)
        return [AMPA_PART, netcon_AMPA], [NMDA_PART, netcon_NMDA]

    def destroy(self):
        pass

res = OPEN_RES()
print(res.get_model())
print(res.res_pos)