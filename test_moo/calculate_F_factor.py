import signal
from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt
from open_pickle import read_from_pickle
from add_figure import add_figure
import pickle
from glob import glob
h.load_file("import3d.hoc")
from math import pi

class Cell: pass
def mkcell(fname):
    #def to read ACS file
    h('objref cell, tobj')
    loader = h.Import3d_GUI(None)
    loader.box.unmap()
    loader.readfile(fname)
    c = Cell()
    loader.instantiate(c)
    return c

def instantiate_swc(filename):
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

######################################################
# build the model
######################################################

fname = "05_08_A_01062017_Splice_shrink_FINISHED_LABEL_Bluecell_spinec91.ASC"
def calculate_F_factor(cell_name,r_head,spine_neck_L,spine_neck_diam,spine_density=1.08,resize_dend=1, shrinkeage=1):
    file_name=glob(cell_name+'*.ASC')[0]
    cell = mkcell(file_name)
    # cell=instantiate_swc('/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/try1.swc')
    for sec in cell.axon:
       h.delete_section(sec=sec)
    for sec in h.allsec():
        sec.diam*=resize_dend
        sec.L*=shrinkeage
    dend_len=np.sum([sec.L for sec in cell.dend])
    spine_in_Micron_density = spine_density  # 12/10 #[N_spine/micrometer] number of spines in micrometer on the dendrite
    head_area=4*pi*r_head**2
    neck_area=2*pi*(spine_neck_diam/2)*spine_neck_L
    spine_area=neck_area+head_area
    spines_area=spine_area*dend_len*spine_in_Micron_density
    dends_area=np.sum([seg.area() for sec in cell.dend for seg in sec]) #* (1.0/0.7)
    F_factor=(spines_area+dends_area)/dends_area
    return F_factor

