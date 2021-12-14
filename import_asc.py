from neuron import h, gui
from matplotlib import pyplot as plt
import numpy as np
from numpy import random as rand
import datetime as dame
import sys
import os.path

h.load_file('import3d.hoc')

arquivo = '/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/171101HuSHS2C1IN0toIN1__postsynaptic_reconstruction_with_putative_synapses.ASC'
arquivo='/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/05_08_A_01062017_Splice_shrink_FINISHED_LABEL_Bluecell_spinec91.ASC'
################~Importing morphology~########################################



def load(arquivo, cell, use_axon=True, xshift=0, yshift=0, zshift=0):
    """load an SWC from filename and instantiate inside cell"""

    name_form = {1: 'soma[%d]', 2: 'axon[%d]', 3: 'dend[%d]', 4: 'apic[%d]'}

    # a helper library, included with NEURON
    h.load_file('stdlib.hoc')
    h.load_file('import3d.hoc')
    h.load_file('stdrun.hoc')

    # load the data. Use Import3d_SWC_read for swc, Import3d_Neurolucida3 for
    # Neurolucida V3, Import3d_MorphML for MorphML (level 1 of NeuroML), or
    # Import3d_Eutectic_read for Eutectic.
    morph = h.Import3d_Neurolucida3()
    morph.input(arquivo)

    # easiest to instantiate by passing the loaded morphology to the Import3d_GUI
    # tool; with a second argument of 0, it won't display the GUI, but it will allow
    # use of the GUI's features
    i3d = h.Import3d_GUI(morph, 0)

    # get a list of the swc section objects
    asc_secs = i3d.asc.sections
    asc_secs = [asc_secs.object(i) for i in np.arange(int(asc_secs.count()))]

    # initialize the lists of sections
    cell.soma, cell.apic, cell.dend, cell.axon = [], [], [], []
    sec_list = {1: cell.soma, 2: cell.axon, 3: cell.dend, 4: cell.apic}

    # name and create the sections
    real_secs = {}
    for asc_sec in asc_secs:
        cell_part = int(asc_sec.type)

        # skip everything else if it's an axon and we're not supposed to
        # use it... or if is_subsidiary
        if (not (use_axon) and cell_part == 2) or asc_sec.is_subsidiary:
            continue

        # figure out the name of the new section
        if cell_part not in name_form:
            raise Exception('unsupported point type')
        name = name_form[cell_part] % len(sec_list[cell_part])

        # create the section
        sec = h.Section(cell=cell, name=name)

        # connect to parent, if any
        if asc_sec.parentsec is not None:
            sec.connect(real_secs[asc_sec.parentsec.hname()](asc_sec.parentx))

        # define shape
        if asc_sec.first == 1:
            h.pt3dstyle(1, asc_sec.raw.getval(0, 0), asc_sec.raw.getval(1, 0),
                        asc_sec.raw.getval(2, 0), sec=sec)

        j = asc_sec.first
        xx, yy, zz = [asc_sec.raw.getrow(i).c(j) for i in np.arange(3)]
        dd = asc_sec.d.c(j)
        if asc_sec.iscontour_:
            # never happens in SWC files, but can happen in other formats supported
            # by NEURON's Import3D GUI
            raise Exception('Unsupported section style: contour')

        if dd.size() == 1:
            # single point soma; treat as sphere
            x, y, z, d = [dim.x[0] for dim in [xx, yy, zz, dd]]
            for xprime in [x - d / 2., x, x + d / 2.]:
                h.pt3dadd(xprime + xshift, y + yshift, z + zshift, d, sec=sec)
        else:
            for x, y, z, d in zip(xx, yy, zz, dd):
                h.pt3dadd(x + xshift, y + yshift, z + zshift, d, sec=sec)

        # store the section in the appropriate list in the cell and lookup table               
        sec_list[cell_part].append(sec)
        real_secs[asc_sec.hname()] = sec

    cell.all = cell.soma + cell.apic + cell.dend + cell.axon


def main(filename=arquivo):
    """demo test program"""

    class Cell:
        def __str__(self):
            return 'neuron'

    cell = Cell()
    load(arquivo, cell)
    return cell


if __name__ == '__main__':
    cell = main()

###############~Visualize imported neuron~####################################
shape_window = h.PlotShape()
shape_window.exec_menu('Show Diam')

input("Press <enter> to continue")


