import signal
from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from add_figure import add_figure
from glob import glob
import os
from find_apic import find_apic
from find_synaptic_loc import syn_dis_from_soma, synaptic_loc

h.load_file("import3d.hoc")
h.load_file("nrngui.hoc")
h.load_file('stdlib.hoc')
h.load_file("stdgui.hoc")


def SIGSEGV_signal_arises(signalNum, stack):
    print(f"{signalNum} : SIGSEGV arises")
    # Your code

signal.signal(signal.SIGSEGV, SIGSEGV_signal_arises)

class Cell() :pass
def mkcell(fname):
    #def to read ACS file
    # morph_reader = h.Import3d_Neurolucida3()
    # morph_reader.input(fname)
    # i3d = h.Import3d_GUI(morph_reader, 0)
    # i3d.instantiate(morph_reader)

    h('objref cell, tobj')
    loader = h.Import3d_GUI(None)
    loader.box.unmap()
    loader.readfile(fname)
    c = Cell()
    loader.instantiate(c)
    return c


def track_one(terminal):
    h.distance(0, 0.5, sec=soma)
    sec=terminal
    dis=[]
    diam=[]
    while sec !=soma:
        if sec.parentseg()==None:
            print("section", sec,"isn't connect",flush=True)
            return [None],[None]
            break

        dis.append(h.distance(sec.parentseg()))
        sec_ref=h.SectionRef(sec=sec)
        diam.append(sec.diam)
        sec=sec_ref.parent
    return np.array(dis),np.array(diam)
######################################################
# build the model
######################################################

yoni_models=glob('*.ASC')


for file in [yoni_models[0]]:
    print('file name: ',file,flush=True)
    cell=[]
    cell=mkcell(file)
    a=[]
    for i in cell.soma:a.append(i)
    if len(a)>1: continue
    print (cell)
    sp = h.PlotShape()
    sp.show(0)
    if file=='05_08_A_01062017_Splice_shrink_FINISHED_LABEL_Bluecell_spinec91.ASC':
        for sec in cell.axon:
            h.delete_section(sec=sec)
        apics = find_apic(cell)
    soma= cell.soma[0]
    # insert pas to all other section
    for sec in tqdm(h.allsec()):
        sec.insert('pas') # insert passive property
        sec.nseg = int(sec.L/10)+1  #decide that the number of segment will be 21 with the same distances


    clamp = h.IClamp(soma(0.5)) # insert clamp(constant potentientiol) at the soma's center
    clamp.amp = +310/1000  ## [nA] supopsed to be 0.05nA
    clamp.dur = 1.1/1000 #[ms]
    clamp.delay = 0.1/1000 #[ms]

    imp = h.Impedance(sec=soma)
    imp.loc(soma(0.5))
    imp.compute(0)
    print('the impadence is',imp.input(0))

    #print the dendrite diameter:
    soma_ref=h.SectionRef(sec=soma)
    print("the soma's childrens diameter is:")
    for i in range(soma_ref.nchild()):
        print(soma_ref.child[i](0).diam)
    length=0
    for dend in h.allsec():
        if dend==soma: continue
        length+=dend.L
    print("total dendritic length is ",length)
    #track from the terminals to the soma

    terminals = []
    for sec in h.allsec():
        if len(sec.children())==0:
            terminals.append(sec)
    plt.close()
    name=file[file.rfind('/')+1:-4]
    add_figure('diam-dis relation \n'+name,'distance from soma','diameter')
    i=0
    for terminal in tqdm(terminals):
        if terminal in cell.dend:
            color='blue'
        else:
            color='red'
        i+=1
        if file == '05_08_A_01062017_Splice_shrink_FINISHED_LABEL_Bluecell_spinec91.ASC':
            if terminal in apics:
                color='red'
        dis,diam=track_one(terminal)
        if dis[0]!=None:
            plt.plot(dis,diam,color,alpha=0.5)
    if file == '05_08_A_01062017_Splice_shrink_FINISHED_LABEL_Bluecell_spinec91.ASC':
        syn_poses = [(-5.56, -325.88, -451.42)]
        synapses_locations = synaptic_loc('05_08_A_01062017', syn_poses)
        syn_dis=syn_dis_from_soma(cell, synapses_locations)
        for syn in syn_dis:
            syn_name = syn[0].name()
            syn_name = syn_name[syn_name.rfind('.')+1:]
            plt.plot(syn[1],syn[0].diam,'*')
            plt.annotate('syn: '+syn_name+'\ndis: '+str(round(syn[1],2)), xy=(syn[1],syn[0].diam), xycoords='data',
                    xytext=(0.8, 0.95), textcoords='axes fraction',
                    arrowprops=dict(facecolor='black', shrink=0.005),
                    horizontalalignment='right', verticalalignment='top',
                    )
        soma_children_diam = [child.diam for child in soma.children()]
    try: os.mkdir('dis_diam')
    except FileExistsError:pass
    plt.savefig('dis_diam/' + name + '.png')
    for sec in h.allsec():
        h.delete_section(sec=sec)
    sp.show(0)
    a=1