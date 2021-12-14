import signal
from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from add_figure import add_figure
from glob import glob
import os


h.load_file("import3d.hoc")
h.load_file("nrngui.hoc")
h.load_file('stdlib.hoc')
h.load_file("stdgui.hoc")
# h.loadfile("stdrun.hoc")

def SIGSEGV_signal_arises(signalNum, stack):
    print(f"{signalNum} : SIGSEGV arises")
    # Your code

signal.signal(signal.SIGSEGV, SIGSEGV_signal_arises)

class Cell() :pass
def mkcell(fname):
    #def to read ACS file
    cell = h.Import3d_SWC_read()
    cell.input(fname)
    i3d = h.Import3d_GUI(cell, 0)
    i3d.instantiate(None)
    return cell

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
# rat_previous=glob('../neuromorpho_data/previous/*/*/*.swc')
# rat_svoboda=glob('../neuromorpho_data/svoboda/*/*.swc')
# mouse=glob('../neuromorpho_data/de paola/*/*.swc')
allen_model=glob('allen_model/*.swc')
one_cell=['allen_model/Oxtr-T2A-Cre_Ai14-314784.06.01.01_656643532_m.swc']

# try: os.mkdir('mouse')
# except FileExistsError: pass
# try: os.mkdir('rat_previous')
# except FileExistsError: pass
# try: os.mkdir('rat_svoboda')
# except FileExistsError: pass
for file in allen_model:
    print('file name: ',file,flush=True)
    cell=mkcell(file)
    a=[]
    for i in h.soma:a.append(i)
    if len(a)>1: continue
    print (cell)
    sp = h.PlotShape()
    sp.show(0)  # show diameters
    for sec in h.axon:
        h.delete_section(sec=sec)

    soma= h.soma[0]
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
    soma_ref=h.SectionRef(sec=h.soma[0])
    print("the soma's childrens diameter is:")
    for i in range(soma_ref.nchild()):
        print(soma_ref.child[i](0).diam)
    length=0
    for dend in h.allsec():
        if dend==h.soma[0]: continue
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
        i+=1
        dis,diam=track_one(terminal)
        if dis[0]!=None:
            plt.plot(dis,diam,alpha=0.5)

    # try: os.mkdir('allen_model/data')
    # except FileExistsError:pass
    # plt.savefig('allen_model/data/'+name+'.png')
    try: os.mkdir(file[:file.rfind('/')]+'/dis_diam')
    except FileExistsError:pass
    plt.savefig(file[:file.rfind('/')]+'/dis_diam/' + name + '.png')
    a=1