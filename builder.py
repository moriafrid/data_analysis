from neuron import h, gui
import signal
class Cell:
    pass

class Builder:
    @staticmethod
    def initiate_cell(fname):
        signal.signal(signal.SIGSEGV, Builder.SIGSEGV_signal_arises)
        h.load_file("import3d.hoc")
        h.load_file("nrngui.hoc")
        h.load_file('stdlib.hoc')
        h.load_file("stdgui.hoc")
        # h.loadfile("stdrun.hoc")
        cell = Builder.mkcell(fname)
        sp = h.PlotShape()
        sp.show(0)  # show diameters
        ## delete all the axons
        for sec in cell.axon:
            h.delete_section(sec=sec)
        for sec in h.allsec():
            sec.insert('pas') # insert passive property
            sec.nseg = int(sec.L/10)+1  #decide that the number of segment will be 21 with the same distances
        return cell

    @staticmethod
    def mkcell(fname):
        #def to read ACS file
        h('objref cell, tobj')
        loader = h.Import3d_GUI(None)
        loader.box.unmap()
        loader.readfile(fname)
        c = Cell()
        loader.instantiate(c)
        return c

    @staticmethod
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

    @staticmethod
    def SIGSEGV_signal_arises(signalNum, stack):
        print(f"{signalNum} : SIGSEGV arises")
        # Your code

    @staticmethod
    def instantiate_unused(filename):
        c = Cell()
        morph_reader = h.Import3d_Neurolucida3()
        morph_reader.input(
            '/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/05_08_A_01062017_Splice_shrink_FINISHED_LABEL_Bluecell_spinec91.ASC')
        i3d = h.Import3d_GUI(morph_reader, 0)
        i3d.insbuilder.initiete_celltantiate(c)

