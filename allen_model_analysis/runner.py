# from pynwb import NWBHDF5IO
# io = NWBHDF5IO('/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/allen_model_analysis/neuronal_model/486262297.nwb', 'r')
# nwbfile_in = io.read()
from allensdk.model.biophys_sim.neuron.hoc_utils import HocUtils
import allensdk.model.biophysical.runner
from allensdk.core.dat_utilities import DatUtilities
from allensdk.core.nwb_data_set import NwbDataSet
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os


class Utils(HocUtils):
    def __init__(self, description):
        super(Utils, self).__init__(description)

    def create_utils(description, model_type=None):
        import logging
        PERISOMATIC_TYPE = "Biophysical - perisomatic"
        ALL_ACTIVE_TYPE = "Biophysical - all active"
        ''' Factory method to create a Utils subclass.

        Parameters
        ----------
        description : Config instance
            used to initialize Utils subclass

        model_type : string
            Must be one of [PERISOMATIC_TYPE, ALL_ACTIVE_TYPE].  If none, defaults to PERISOMATIC_TYPE

        Returns
        -------
        Utils instance
        '''

        if model_type is None:
            try:
                model_type = description.data['biophys'][0]['model_type']
            except KeyError as e:
                logging.error("Could not infer model type from description")

        if model_type == PERISOMATIC_TYPE:
            return Utils(description)
        elif model_type == ALL_ACTIVE_TYPE:
            return AllActiveUtils(description)


class AllActiveUtils(Utils):
    def generate_morphology(self, morph_filename):
        '''Load a neurolucida or swc-format cell morphology file.

        Parameters
        ----------
        morph_filename : string
            Path to morphology.
        '''
        morph_basename = os.path.basename(morph_filename)
        morph_extension = morph_basename.split('.')[-1]
        if morph_extension.lower() == 'swc':
            morph = self.h.Import3d_SWC_read()
        elif morph_extension.lower() == 'asc':
            morph = self.h.Import3d_Neurolucida3()
        else:
            raise Exception("Unknown filetype: %s" % morph_extension)

        morph.input(morph_filename)
        imprt = self.h.Import3d_GUI(morph, 0)

        self.h("objref this")
        imprt.instantiate(self.h.this)

        for sec in self.h.allsec():
            sec.nseg = 1 + 2 * int(sec.L / 40.0)

        self.h("soma[0] area(0.5)")
        axon_diams = [self.h.axon[0].diam, self.h.axon[0].diam]
        self.h.distance(sec=self.h.soma[0])
        for sec in self.h.allsec():
            if sec.name()[:4] == "axon":
                if self.h.distance(0.5, sec=sec) > 60:
                    axon_diams[1] = sec.diam
                    break
        for sec in self.h.allsec():
            if sec.name()[:4] == "axon":
                self.h.delete_section(sec=sec)
        self.h('create axon[2]')
        for index, sec in enumerate(self.h.axon):
            sec.L = 30
            sec.diam = axon_diams[index]

        for sec in self.h.allsec():
            sec.nseg = 1 + 2 * int(sec.L / 40.0)

        self.h.axon[0].connect(self.h.soma[0], 1.0, 0.0)
        self.h.axon[1].connect(self.h.axon[0], 1.0, 0.0)

        # make sure diam reflects 3d points
        self.h.area(.5, sec=self.h.soma[0])

def load_cell_parameters(self):
    '''Configure a neuron after the cell morphology has been loaded.'''
    passive = self.description.data['passive'][0]
    genome = self.description.data['genome']
    conditions = self.description.data['conditions'][0]
    h = self.h

    h("access soma")

    # Set fixed passive properties
    for sec in h.allsec():
        sec.Ra = passive['ra']
        sec.insert('pas')
        # for seg in sec:
        #     seg.pas.e = passive["e_pas"]

    # for c in passive["cm"]:
    #     h('forsec "' + c["section"] + '" { cm = %g }' % c["cm"])

    # Insert channels and set parameters
    for p in genome:
        section_array = p["section"]
        mechanism = p["mechanism"]
        param_name = p["name"]
        param_value = float(p["value"])
        if section_array == "glob":  # global parameter
            h(p["name"] + " = %g " % p["value"])
        else:
            if hasattr(h, section_array):
                if mechanism != "":
                    print('Adding mechanism %s to %s'
                          % (mechanism, section_array))
                    for section in getattr(h, section_array):
                        if self.h.ismembrane(str(mechanism),
                                             sec=section) != 1:
                            section.insert(mechanism)

                print('Setting %s to %.6g in %s'
                      % (param_name, param_value, section_array))
                for section in getattr(h, section_array):
                    setattr(section, param_name, param_value)

    # Set reversal potentials
    for erev in conditions['erev']:
        erev_section_array = erev["section"]
        ek = float(erev["ek"])
        ena = float(erev["ena"])

        print('Setting ek to %.6g and ena to %.6g in %s'
              % (ek, ena, erev_section_array))

        if hasattr(h, erev_section_array):
            for section in getattr(h, erev_section_array):
                if self.h.ismembrane("k_ion", sec=section) == 1:
                    setattr(section, 'ek', ek)

                if self.h.ismembrane("na_ion", sec=section) == 1:
                    setattr(section, 'ena', ena)
        else:
            print("Warning: can't set erev for %s, "
                  "section array doesn't exist" % erev_section_array)

    self.h.v_init = conditions['v_init']
    self.h.celsius = conditions['celsius']



#description="Biophysical - all active"
description='/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/allen_model_analysis/neuronal_model/manifest.json'
utils = Utils.create_utils(description,model_type="Biophysical - all active")
h = utils.h

# configure model
manifest = description.manifest
morphology_path = description.manifest.get_path('Cux2-CreERT2_Ai14-207761.04.02.01_506798042_m.swc')
utils.generate_morphology(morphology_path.encode('ascii', 'ignore'))
utils.load_cell_parameters()

# configure stimulus and recording
stimulus_path = description.manifest.get_path('stimulus_path')
nwb_out_path = manifest.get_path("output")
output = NwbDataSet(nwb_out_path)
run_params = description.data['runs'][0]
sweeps = run_params['sweeps']
junction_potential = description.data['fitting'][0]['junction_potential']
mV = 1.0e-3
# run sweeps
for sweep in sweeps:
    utils.setup_iclamp(stimulus_path, sweep=sweep)
    vec = utils.record_values()

    h.finitialize()
    h.run()

    # write to an NWB File
    output_data = (np.array(vec['v']) - junction_potential) * mV
    output.set_sweep(sweep, None, output_data)


nwb_file = '486262297_ephys.nwb'
sweep_number = 52
dat_file = '486262297_ephys_%d.dat' % (sweep_number)

nwb = NwbDataSet(nwb_file)
sweep = nwb.get_sweep(sweep_number)

# read v and t as numpy arrays
v = sweep['response']
dt = 1.0e3 / sweep['sampling_rate']
num_samples = len(v)
t = np.arange(num_samples) * dt

# save as text file
data = np.transpose(np.vstack((t, v)))
with open (dat_file, "w") as f:
    np.savetxt(f, data)

# save image using matplotlib
fig, ax = plt.subplots(nrows=1, ncols=1)
ax.plot(t, v)
ax.set_title("Sweep %s" % (sweep_number))
fig.savefig('out.png')

stimulus_path = description.manifest.get_path('stimulus_path')
nwb_out_path = manifest.get_path("output")
output = NwbDataSet(nwb_out_path)
run_params = description.data['runs'][0]
sweeps = run_params['sweeps']
junction_potential = description.data['fitting'][0]['junction_potential']
mV = 1.0e-3