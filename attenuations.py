from neuron import h, gui
from glob import glob
import numpy as np
import os
import matplotlib.pyplot as plt
import signal
freq=100
resize_diam_by=2
do_resize_dend=True
norm_Rin=False
folder_='/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/'
h.load_file("import3d.hoc")
h.load_file("nrngui.hoc")
h.load_file('stdlib.hoc')
h.load_file("stdgui.hoc")

cell_name='05_08_A_01062017'
class Cell: pass
def mkcell(fname):
    #def to read ACS file
  loader = h.Import3d_GUI(None)
  loader.box.unmap()
  loader.readfile(fname)
  c = Cell()
  loader.instantiate(c)
  return c

def instantiate_swc(filename):
    h.load_file('import3d.hoc')
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
def change_model_pas(cell, CM=1, RA = 250, RM = 20000.0, E_PAS = -77.5, F_factor = {}, SPINE_START=60):
    #input the neuron property    h.dt = 0.1

    h.distance(0,0.5, sec=cell.soma[0]) # it isn't good beacause it change the synapse distance to the soma
    #h.distance(0, sec=soma)
    for sec in cell.all: ##check if we need to insert Ra,cm,g_pas,e_pas to the dendrit or just to the soma
        if isinstance(cell, Cell):
            if sec in cell.axon: continue   #@# cell.axon is not exist in hoc object
        sec.Ra = RA
        sec.cm = CM
        sec.g_pas = 1.0 / RM
        sec.e_pas = E_PAS
    for sec in cell.dend:
        for seg in sec:
            if h.distance(seg) > SPINE_START:
                if type(F_factor)!=float :
                    F_factor = 2.0 # F_factor[sec]
                seg.cm *= F_factor
                seg.g_pas *= F_factor

def plot_records(RM, RA, CM,cell, syn,spine=None,save_name= "lambda"):
    soma = cell.soma[0]
    # syn = cell.dend[82]
    folder_="data/Attenuations/"
    try: os.mkdir(folder_)
    except FileExistsError: pass
    if syn_injection:
        try: os.mkdir(folder_+'syn/')
        except FileExistsError: pass
    elif clamp_injection:
        try:os.mkdir(folder_ + 'clamp_inj/')
        except FileExistsError:pass

    change_model_pas(cell,CM=CM, RA=RA, RM=RM, E_PAS = E_PAS,F_factor=F_factor_result)
    Vvec_soma = h.Vector()
    Vvec_soma.record(soma(0.5)._ref_v)

    Vvec_dend = h.Vector()
    Vvec_dend.record(syn[0](syn[1])._ref_v)
    Vvec_spine = h.Vector()

    Tvec = h.Vector()
    Tvec.record(h._ref_t)

    if clamp_injection:
        Ivec = h.Vector()
        Ivec.record(clamp._ref_i)
        if spine!=None:
            Vvec_spine.record(spine(1)._ref_v)


    if syn_injection:
        if spine!=None:
            Vvec_spine.record(spine(1)._ref_v)
            Ivec=Vvec_spine
        else:
            h.run()
            Ivec=range(len(Tvec))
        # Ivec.record(netcon._ref_i)
    h.run()
    npTvec = np.array(Tvec)
    npIvec = np.array(Ivec)
    npVec_soma = np.array(Vvec_soma)
    npVec_dend = np.array(Vvec_dend)
    if norm_Rin:
        npIvec/=Rin_syn_resize_dend #npIvec*(Rin_syn_0/Rin_syn_resize_dend)
        # npVec_dend=npVec_dend*(Rin_dend_0/Rin_dend_resize_dend)#/=Rin_dend_resize_dend #npVec_dend*(Rin_dend_0/Rin_dend_resize_dend)
        # npVec_soma=npVec_soma*(Rin_soma_0/Rin_soma_resize_dend)#/=Rin_soma_resize_dend#npVec_soma*(Rin_soma_0/Rin_soma_resize_dend)
    from add_figure import add_figure
    figure, axis = plt.subplots(3, 1)
    axis[0].plot(npTvec,npIvec)
    # axis[0].set_title("spine voltage")
    axis[0].set_xlabel('mS')
    axis[0].set_ylabel('ms')

    axis[1].plot(npTvec, npVec_dend)
    axis[1].set_title("dend parent spine voltage")
    axis[1].set_xlabel('mS')
    axis[1].set_ylabel('mV')

    axis[2].plot(npTvec, npVec_soma)
    axis[2].set_title("soma Voltage")
    axis[2].set_xlabel('mS')
    axis[2].set_ylabel('mv')
    if clamp_injection:

        axis[0].set_title("\n current injection of " + str(clamp.amp) + "nA to the syn for " + str(pulse_size) + 'ms')
        figure.tight_layout(pad=1.0)

        if do_resize_dend:
            plt.savefig(folder_ + 'clamp_inj freq_'+str(freq)+'/' + str(pulse_size) + "ms_dend*"+str(resize_diam_by)+'.png')


    elif syn_injection:
        axis[0].set_title("syn weight of " + str(syn_weight) + '\nspine head Volt/Rinput')
        figure.tight_layout(pad=1.0)
        if not norm_Rin:
            axis[0].set_title("syn weight of " + str(syn_weight) + '\nspine head Voltage')
            axis[0].set_ylabel('mv')
        if do_resize_dend:
            plt.savefig(folder_+'syn'++'/'+str(syn_weight)+"_weight_dend*"+str(resize_diam_by)+".png")
        else:
            plt.savefig(folder_+'syn/'+str(syn_weight)+"_weight.png")



def create_spine( icell, sec, pos, number=0, neck_diam=0.25, neck_length=1.35,
                 head_diam=0.944):  # np.sqrt(2.8/(4*np.pi))
    neck = h.Section(name="spineNeck" + str(number))
    head = h.Section(name="spineHead" + str(number))
    neck.L = neck_length
    neck.diam = neck_diam
    head.diam = head_diam
    head.L = head_diam
    head.connect(neck(1))
    neck.connect(sec(pos))
    h("access " + str(neck.hoc_internal_name()))
    icell.all.append(neck)
    # icell.dend.append(neck)

    # if Rneck == "normal_neck":
    #     icell.all.append(neck)
    #     if sec.name().find('dend') > -1:
    #         icell.basal.append(neck)
    #     else:
    #         icell.apical.append(neck)
    h.pop_section()
    h("access " + str(head.hoc_internal_name()))
    icell.all.append(head)
    # icell.dend.append(head)

    # if sec.name().find('dend') > -1:
    #     icell.basal.append(head)
    # else:
    #     icell.apical.append(head)
    # sim.neuron.h.pop_section()
    h.pop_section()

    for sec in [neck, head]:
        sec.insert("pas")
    if not Rneck == "normal_neck":
        neck.g_pas = 1.0 / passive_val[cell_name]["RM"]
        neck.cm= passive_val[cell_name]["CM"]
        neck.Ra=passive_val[cell_name]["RA"]#int(Rneck)
    return icell,[neck, head]

def add_morph(cell, syn,spine_property,number=0):
    # sim.neuron.h.execute('create spineNeck['+str(len(syns))+']', icell)
    # sim.neuron.h.execute('create spineHead['+str(len(syns))+']', icell)
    cell,spine=create_spine(cell, syn[0],syn[1] ,number=0, neck_diam=spine_property['NECK_DIAM'], neck_length=spine_property['NECK_LENGHT'],head_diam=spine_property['HEAD_DIAM'])
    return cell,spine
    # num = syn[0]
        # num = int(num[num.find("[") + 1:num.find("]")])
        # if syn[0].find("dend") > -1:
        #     sec = cell.dend[num]
        # elif syn[0].find("apic") > -1:
        #     sec = cell.apic[num]
        # else:
        #     sec = cell.soma[0]
        # all.append(create_spine(cell, sec, syn[0],syn[1] ,number=i, neck_diam=spine_property[str(i)]['NECK_DIAM']), neck_length=spine_property[str(i)]['NECK_LENGHT'],
        #                         head_diam=spine_property[str(i)]['HEAD_DIAM'])
    # return all
# cell=instantiate_swc('/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/try1.swc')
fname = glob(cell_name+'*.ASC')[0]
cell =mkcell(fname)
print (cell)
for sec in cell.axon:
   h.delete_section(sec=sec)
soma = cell.soma[0]
from find_synaptic_loc import synaptic_loc
syn_poses={}
syn_poses['05_08_A_01062017']=[(-5.56, -325.88, -451.42)]
syns = synaptic_loc(cell,syn_poses[cell_name],del_axon=False)['place_as_sec']

for sec in h.allsec():
    # if isinstance(cell, Cell):
    # if sec in cell.axon: continue   #@# cell.axon is not exist in hoc object
    sec.insert('pas') # insert passive property
    sec.nseg = int(sec.L/10)+1  #decide that the number of segment will be 21 with the same distances
syn=syns[0]
# if do_resize_dend:
#     imp_0 = h.Impedance(sec=syn[0])
#     imp_0.loc(0.165, sec=syn[0])
#     imp_0.compute(0)  # check if you need at 10 Hz
#     Rin_syn_0 = imp_0.input(syn[1], sec=syn[0])
# # for sec in cell.dend:
# #     sec.diam = sec.diam*resize_diam_by
# if norm_Rin:
#     imp=h.Impedance(sec=syn[0])
#     imp.loc(syn[1], sec=syn[0])
#     imp.compute(0) #check if you need at 10 Hz
#     Rin_syn_resize_dend = imp.input(syn[1], sec=syn[0])
#@# I need to think on a way to do it more clever

# for i,syn in enumerate(syns):
#     spine_parameters[str(i)]={'location':[syn],'NECK_LENGHT':0.782,'spine_head_vol':}
passive_val = {'05_08_A_01062017': {'CM': 1.88, 'RA': 95.7, 'RM': 12371}}
from math import pi
from calculate_F_factor import calculate_F_factor
if cell_name=='05_08_A_01062017':
    NECK_LENGHT=0.782
    spine_head_vol=0.139
    R_head=(spine_head_vol/(4*pi/3))**(1/3) #0.32
    spine_neck_diam = 0.164
    F_factor_result=calculate_F_factor(cell, spine_head_vol,spine_neck_diam, NECK_LENGHT)
    HEAD_DIAM=2*R_head#0.64
    Rneck = passive_val[cell_name]["RA"]
spines_property={'0':{'NECK_LENGHT':NECK_LENGHT,'NECK_DIAM':spine_neck_diam,'HEAD_DIAM':HEAD_DIAM}}
sp = h.PlotShape()
sp.show(0)  # show diameters
sp.color(2, sec=cell.dend[82] )
pulse_size=1000
clamp_injection=True
syn_injection=False
syn_on_dend=False
put_syn_on_spine_head=False
# h.celsius = 36

if clamp_injection:
    clamp = h.IClamp(cell.dend[82](0.992)) # insert clamp(constant potentientiol) at the soma's center
    clamp.amp = 0.05 #nA
    clamp.delay = 200 #ms
    clamp.dur = pulse_size  # ms
    h.tstop = 500

elif syn_injection: #@# replace in my specipic synapse
    # create excitatory synapse at knowen head spine
    spines=[]
    syn_shape=[]
    syn=syns[0]
    if put_syn_on_spine_head:
        for i,syn in enumerate(syns):
            cell,spine=add_morph(cell, syn ,spines_property[str(i)])
            spines.append(spine)
            syn_shape = h.Exp2Syn(spine[0](1))
    else:
        syn_shape=h.Exp2Syn(syn[0](syn[1]))
    syn_shape.tau1 = 0.3
    syn_shape.tau2 = 1.8
    syn_shape.e = 0
    # syn_shape.i=0
    # syn_shape.onset=200
    # syn_shape.imax=0.5
    syn_netstim = h.NetStim()  # the location of the NetStim does not matter
    syn_netstim.number = 1
    syn_netstim.start = 200
    syn_netstim.noise = 0
    netcon = h.NetCon(syn_netstim, syn_shape)
    syn_weight = 0.002
    netcon.weight[0] = syn_weight  # activate the on_path synapse
    # create a NetStim that will activate the synapses at t=1000
    h.tstop = 1000

    ####### if there is more then 1 synapse
    # syn_shape=[]
    # for i, syn in enumerate(syns):
    #     syn_shape.append(h.Exp2Syn(syn[0](syn[1])))
    #     syn_shape[i].tau1 = 0.3
    #     syn_shape[i].tau2 = 1.8
    #     syn_shape[i].e = 0
    #
    #     # syn_shape.i=0
    #     # syn_shape.onset=200
    #     # syn_shape.imax=0.5
    #     syn_netstim = h.NetStim()  # the location of the NetStim does not matter
    #     syn_netstim.number = 1
    #     syn_netstim.start = 200
    #     syn_netstim.noise = 0
    #     netcon = h.NetCon(syn_netstim, syn_shape[i])
    #     syn_weight = 0.003
    #     netcon.weight[i] = syn_weight  # activate the on_path synapse
    #     # create a NetStim that will activate the synapses  at t=1000
    #
    # # syn_weight=0.003
    # # netcon.weight = syn_weight  # activate the on_path synapse
    # # h.tstop = 1000
if do_resize_dend:
    imp = h.Impedance(sec=spine[1])
    imp.loc(1, sec=spine[1])
    imp.compute(freq)  # check if you need at 10 Hz
    Rin_syn_0 = imp.input(1, spine[1])

    imp = h.Impedance(sec=syn[0])
    imp.loc(syn[1], sec=syn[0])
    imp.compute(freq)  # check if you need at 10 Hz
    Rin_dend_0 = imp.input(syn[1], sec=syn[0])

    imp = h.Impedance(sec=soma)
    imp.loc(syn[1], sec=soma)
    imp.compute(freq)  # check if you need at 10 Hz
    Rin_soma_0 = imp.input(0.5, sec=soma)

    for sec in cell.dend:
        sec.diam = sec.diam*resize_diam_by
if norm_Rin:
    syn=syns[0]
    imp=h.Impedance(sec=spine[1])
    imp.loc(1, sec=spine[1])
    imp.compute(freq) #check if you need at 10 Hz
    Rin_syn_resize_dend = imp.input(1, sec=spine[1])

    imp=h.Impedance(sec=syn[0])
    imp.loc(syn[1], sec=syn[0])
    imp.compute(freq) #check if you need at 10 Hz
    Rin_dend_resize_dend = imp.input(syn[1], sec=syn[0])

    imp = h.Impedance(sec=soma)
    imp.loc(syn[1], sec=soma)
    imp.compute(freq)  # check if you need at 10 Hz
    Rin_soma_resize_dend = imp.input(0.5, sec=soma)
E_PAS=-77.5

h.v_init=E_PAS
h.dt = 0.1 #=hz
h.steps_per_ms = h.dt

RM_const = 60000.0
RA_const = 250.0
CM_const = 1.0

CM=2/2
RM=5684*2 #20000 #5684*2
RA=100
if put_syn_on_spine_head:
    plot_records(RM, RA, CM,cell,syns[0], spine=spine[1],save_name= "lambda")
else:
    plot_records(RM, RA, CM,cell,syns[0],save_name= "lambda")


a=1