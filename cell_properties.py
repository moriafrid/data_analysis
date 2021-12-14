import signal
from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from add_figure import add_figure

path='/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/data/short_pulse_mean.p'
folder_='/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/'
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
  loader = h.Import3d_GUI(None)
  loader.box.unmap()
  loader.readfile(fname)
  c = Cell()
  loader.instantiate(c)
  return c

######################################################
# build the model
######################################################

fname = "05_08_A_01062017_Splice_shrink_FINISHED_LABEL_Bluecell_spinec91.ASC"
cell=mkcell(fname)
print (cell)
sp = h.PlotShape()
sp.show(0)  # show diameters

for sec in cell.axon:
   h.delete_section(sec=sec)

soma= cell.soma[0]
# insert pas to all other section
for sec in tqdm(h.allsec()):
    sec.insert('pas') # insert passive property
    sec.nseg = int(sec.L/10)+1  #decide that the number of segment will be 21 with the same distances


clamp = h.IClamp(soma(0.5)) # insert clamp(constant potentientiol) at the soma's center
clamp.amp = -0.05  ## supopsed to be 0.05nA
clamp.dur = 200
clamp.delay = 190

#short_pulse, hz, rest = read_from_pickle(path, hz=True ,rest=True)
# V = short_pulse[0]
# T = short_pulse[1]
# E_PAS = rest
# V+=E_PAS
#
# h.tstop = T[-1]
# h.v_init=E_PAS
# h.dt = 0.1
# h.steps_per_ms = h.dt

imp = h.Impedance(sec=soma)
imp.loc(soma(0.5))
imp.compute(0)
imp.input(0)
imp.compute(0)
print('the impadence is',imp.input(0))

#print the dendrite diameter:
soma_ref=h.SectionRef(sec=cell.soma[0])
print("the soma's childrens diameter is:")
for i in range(soma_ref.nchild()):
    print(soma_ref.child[i](0).diam)
length=0
for dend in cell.dend:
    length+=dend.L
print("total dendritic length is ",length)
#track from the terminals to the soma
def track_one(terminal):
    h.distance(0, 0.5, sec=soma)
    sec=terminal
    dis=[]
    diam=[]
    while sec !=soma:
        dis.append(h.distance(sec.parentseg()))
        sec_ref=h.SectionRef(sec=sec)
        diam.append(sec.diam)
        sec=sec_ref.parent
    return np.array(dis),np.array(diam)
terminals = []
for sec in cell.dend:
    if len(sec.children())==0:
        terminals.append(sec)
plt.close()
add_figure('diam-dis relation along dendrites with diffrent collors','distance from soma','diameter')
i=0
for terminal in terminals[:-14:2]:
    i+=1
    dis,diam=track_one(terminal)
    plt.plot(dis,diam,alpha=0.5)
plt.savefig('data/diam-dis.png')
plt.savefig('data/diam-dis.pdf')
plt.show()
a=1