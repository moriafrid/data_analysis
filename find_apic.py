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


def SIGSEGV_signal_arises(signalNum, stack):
    print(f"{signalNum} : SIGSEGV arises")
    # Your code

signal.signal(signal.SIGSEGV, SIGSEGV_signal_arises)

class Cell() :pass
def mkcell(fname):
    h('objref cell, tobj')
    loader = h.Import3d_GUI(None)
    loader.box.unmap()
    loader.readfile(fname)
    c = Cell()
    loader.instantiate(c)
    return c
def run(apics,last_apic):
    for child in last_apic.children():
        apics.append(child)
        run(apics,child)
def find_apic(cell,del_axon=False):
    print("unsure that there isn't cell.axon that send to find_apic.py or del_axon=False")
    last_dend_diam=0
    for dend in cell.soma[0].children():
        if not del_axon:
            if dend in cell.axon: continue
        if dend.diam>last_dend_diam:
            start_apic=dend
            last_dend_diam=dend.diam
    apics=[start_apic]
    run(apics,start_apic)
    return apics







if __name__=='__main__':
    fname='05_08_A_01062017_Splice_shrink_FINISHED_LABEL_Bluecell_spinec91.ASC'
    cell=mkcell(fname)
    apics=find_apic(cell)
    a=1
