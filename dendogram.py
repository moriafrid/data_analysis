import numpy as np
from neuron import h, gui
import sys,os
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from glob import glob
import signal
from find_apic import find_apic
from find_synaptic_loc import synaptic_loc

h.load_file("import3d.hoc")
h.load_file("nrngui.hoc")
h.load_file('stdlib.hoc')
h.load_file("stdgui.hoc")
colors_dict = {"soma":"k",
               "apical": "blue",
               "oblique":"cyan",
               "trunk":"purple",
               "basal": "r",
               "axon": "green",
               "else": "gold",
               "synapse": "grey"}

def SIGSEGV_signal_arises(signalNum, stack):
    print(f"{signalNum} : SIGSEGV arises")
    # Your code

signal.signal(signal.SIGSEGV, SIGSEGV_signal_arises)

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

def get_segment_length_lamda(seg):
    """
	return the segment  e_length
	:param seg_len:
	:param RM:
	:param RA:
	:return:
	"""
    sec = seg.sec
    seg_len = sec.L/sec.nseg #micro meter
    d = seg.diam #micro meter
    R_total = 1.0 / seg.g_pas #Rm[cm^2*oum] sce.Ra[cm*oum]
    lamda = np.sqrt((R_total / sec.Ra) * (d / 10000.0) / 4.0) #micro meter
    # return (float(seg_len) / 10000.0) / lamda
    return lamda



def add_sec(self, sec):
    """
    electric dendogram
    :param sec:
    :return:
    """
    sec_length = 0
    for seg in sec:
        sec_length += get_segment_length_lamda(seg)
    parent = h.SectionRef(sec=sec).parent
    self.tree_dendogram_dist[sec] = self.tree_dendogram_dist[parent] + sec_length


def add_sec2(self, sec):
    """
    morpho dendogram
    :param sec:
    :return:
    """
    h.distance(0, 0.5, sec=self.cell.soma[0])
    self.tree_dendogram_dist[sec] = h.distance(1, sec=sec)

def get_spine_area():
    neck_length=0.78
    neck_diam = 1.64
    head_volume = 0.14
    head_r = (head_volume*3/4/np.pi)**(1/3)
    head_area = 4*np.pi*head_r**3
    neck_area = np.pi * neck_diam * neck_length
    return head_area +neck_area

# def instantiate_swc(filename):
#     h('objref cell, tobj')
#     h.load_file('allen_model.hoc')
#     h.execute('cell = new allen_model()')
#     h.load_file(filename)
#     nl = h.Import3d_SWC_read()
#     nl.quiet = 1
#     nl.input(filename)
#     i3d = h.Import3d_GUI(nl, 0)
#     i3d.instantiate(h.cell)
#     return h.cell

def change_model_pas(cell, CM=1, RA = 250, RM = 20000.0, E_PAS = -70.0, F_factor = {}):
    #input the neuron property    h.dt = 0.1

    h.distance(0,0.5, sec=cell.soma[0]) # it isn't good beacause it change the synapse distance to the soma
    for sec in h.allsec(): ##check if we need to insert Ra,cm,g_pas,e_pas to the dendrit or just to the soma
        sec.Ra = RA
        sec.cm = CM
        sec.g_pas = 1.0 / RM
        sec.e_pas = E_PAS
    for sec in cell.dend:
        for seg in sec: #count the number of segment and calclate g_factor and total dend distance,
            # how many segment have diffrent space larger then SPINE_START that decided
            if h.distance(seg) > 60:
                F_factor = F_factor#2.03 # F_factor[sec]
                seg.cm *= F_factor
                seg.g_pas *= F_factor
    return cell

class Dendogram():
    def __init__(self,
                 name,
                 morph_path,
                 length_function,
                 dots_loc=[],
                 color_dict = colors_dict,
                 diam_factor=None,
                 del_axon=True,
                 apic=[],
                 F_factor=2.03):
        self.name=name
        self.colors_dict = color_dict
        self.cell = mkcell(morph_path)
        if del_axon:
            for sec in self.cell.axon:
                h.delete_section(sec=sec)
        for sec in h.allsec():
            sec.insert('pas')
            sec.nseg = max(int(sec.L), 1)

        self.cell=change_model_pas(self.cell, CM=1.88, RA = 95.7, RM = 12371, E_PAS = -77.0,F_factor=2.03)
        self.morph_path=morph_path
        self.tree_dendogram_dist = dict()
        self.tree_dendogram_dist[self.cell.soma[0]] = 0
        self.add_sec = length_function
        self.diam_factor=diam_factor
        self.apic= find_apic(self.cell,del_axon=del_axon)
        self.syn = synaptic_loc(self.cell,dots_loc)['place_name']
        self.dots_loc = np.array([[syn[0],syn[1]] for syn in self.syn])
        # self.dots_loc = np.array(dots_loc)
        # self.dots_loc = np.array([[eval('self.cell.'+syn_[0], {'self':self}, {}),syn_[1]] for syn_ in self.syn])

    def cumpute_distances(self, base_sec):
        for sec in h.SectionRef(sec=base_sec).child:
            self.add_sec(self, sec)
            self.cumpute_distances(sec)

    def add_syn(self, sec, seg):
        h.distance(0, 0.5, sec=self.cell.soma)
        self.tree_dendogram_dist[sec] = h.distance(seg, sec=sec)

    def get_color(self, sec):
        if sec in self.cell.dend:
            if sec in self.apic:
                return self.colors_dict["apical"]
            else:
                return self.colors_dict["basal"]

        # elif sec in self.cell.apic:
        #     return self.colors_dict["apical"]

        elif sec in self.cell.axon:
            return self.colors_dict["axon"]
        elif sec in self.cell.soma:
            return self.colors_dict["soma"]
        else:
            return self.colors_dict["else"]

    def plot_synapse(self, sec_start, sec_end, pos, x_axis):
        syn_dis=sec_start + abs(sec_end - sec_start) * float(pos)
        plt.scatter(x_axis, syn_dis, color=colors_dict["synapse"])
        plt.annotate('syn distance\n' + str(round(syn_dis,3)),fontsize=10, xy=(x_axis, syn_dis), xycoords='data',
                     xytext=(0.3, 0.5), textcoords='axes fraction',
                     arrowprops=dict(facecolor='black', shrink=0.005,lw=0.0005),
                     horizontalalignment='center', verticalalignment='top',
                     )
        return syn_dis
    def plot_func(self, sec, x_pos, color):
        parent = h.SectionRef(sec=sec).parent

        if sec in self.done_section:
            raise BaseException("problem with morph")
        else:
            self.done_section.add(sec)
        sec_name = sec.name()
        sec_name = sec_name[sec_name.find(".") + 1:]
        sec_name = sec_name[sec_name.find(".") + 1:]

        if h.SectionRef(sec=sec).nchild() == 0:
            plt.plot([x_pos, x_pos], [self.tree_dendogram_dist[parent], self.tree_dendogram_dist[sec]], color=self.get_color(sec),
                     linewidth=1 if self.diam_factor is None else sec.diam*self.diam_factor)
            for synapse in self.dots_loc:
                sec_n=synapse[0]
                loc=synapse[1]
                if sec_name == sec_n:
                    dendogram_syn_distance =self.plot_synapse(self.tree_dendogram_dist[parent], self.tree_dendogram_dist[sec], loc, x_pos)
                    print('synapse dendogram len is ',dendogram_syn_distance)
            return x_pos + 1.0, x_pos

        elif h.SectionRef(sec=sec).nchild() == 1:
            for synapse in self.dots_loc:
                sec_n=synapse[0]
                loc=synapse[1]
                if sec_name == sec_n:
                    self.plot_synapse(self.tree_dendogram_dist[parent], self.tree_dendogram_dist[sec], loc, x_pos)
            # x_pos+=1
            x_pos, start_pos = self.plot_func(h.SectionRef(sec=sec).child[0], x_pos, color)
            plt.plot([start_pos, start_pos], [self.tree_dendogram_dist[parent], self.tree_dendogram_dist[sec]], color=self.get_color(sec),
                     linewidth=1 if self.diam_factor is None else sec.diam*self.diam_factor)
            return x_pos, start_pos

        x_pos, start_pos = self.plot_func(h.SectionRef(sec=sec).child[0], x_pos, color)
        for i in range(1, int(h.SectionRef(sec=sec).nchild()) - 1, 1):
            x_pos, end_pos = self.plot_func(h.SectionRef(sec=sec).child[i], x_pos, color)

        x_pos, end_pos = self.plot_func(h.SectionRef(sec=sec).child[int(h.SectionRef(sec=sec).nchild()) - 1], x_pos,
                                   color)
        mid_x = start_pos + abs(end_pos - start_pos) / 2.0
        plt.plot([mid_x, mid_x], [self.tree_dendogram_dist[parent], self.tree_dendogram_dist[sec]], color=self.get_color(sec),
                 linewidth=1 if self.diam_factor is None else sec.diam*self.diam_factor)
        plt.plot([start_pos, end_pos], [self.tree_dendogram_dist[sec]] * 2, color=self.get_color(sec), linewidth=1 if self.diam_factor is None else sec.diam*self.diam_factor)

        for sec_n, loc in self.dots_loc:
            if sec_name == sec_n:
                self.plot_synapse(self.tree_dendogram_dist[parent], self.tree_dendogram_dist[sec], loc, mid_x)

        return x_pos, mid_x

    def plot(self, save_folder, max_y=None,title='Dendogram',ylabel='distance from soma'):
        plt.figure(figsize=(10, 10))
        plt.title(title,fontsize=24)
        plt.ylabel(ylabel,fontsize=16)
        x_pos = 0.0
        start_pos=0.0
        self.done_section = set()
        for i in range(0, int(h.SectionRef(sec=self.cell.soma[0]).nchild()), 1):
            sec = h.SectionRef(sec=self.cell.soma[0]).child[i]
            if sec in self.apic:
                x_pos, start_pos = self.plot_func(sec, x_pos, color=self.get_color(sec))
        for i in range(0, int(h.SectionRef(sec=self.cell.soma[0]).nchild()), 1):
            sec = h.SectionRef(sec=self.cell.soma[0]).child[i]
            if sec not in self.apic:
                x_pos, end_pos = self.plot_func(sec, x_pos, color=self.get_color(sec))
        plt.plot([start_pos, end_pos], [0] * 2, color=self.colors_dict["soma"], linewidth=1 if self.diam_factor is None else self.cell.soma[0].diam *self.diam_factor)
        mid_x = start_pos + abs(end_pos - start_pos) / 2.0
        plt.plot([mid_x, mid_x], [-0.01, 0], color=self.colors_dict["soma"], linewidth=1 if self.diam_factor is None else self.cell.soma[0].diam *self.diam_factor)
        plt.xticks([])

        legend_elements = [
            Line2D([0], [0], color=self.colors_dict["soma"], lw=2, label="soma"),
            Line2D([0], [0], color=self.colors_dict["apical"], lw=2, label="apical"),
            Line2D([0], [0], color=self.colors_dict["basal"], lw=2, label="basal"),
            Line2D([0], [0], color=self.colors_dict["trunk"], lw=2, label="trunk"),
            Line2D([0], [0], color=self.colors_dict["oblique"], lw=2, label="oblique"),
            Line2D([0], [0], color=self.colors_dict["synapse"], lw=2, label="synapse")
        ]
        plt.legend(handles=legend_elements, loc="best")
        if max_y is None:
            max_y = plt.ylim()[1]
        plt.ylim([-0.1, max_y])
        plt.savefig(save_folder + self.name)
        plt.savefig(save_folder + self.name+ ".pdf")
        plt.close()
        self.done_section = set()
        return max_y


save_folder = 'E_dendogram/'
save_folder2 = 'M_dendogram/'
try: os.mkdir(save_folder)
except:pass
try: os.mkdir(save_folder2)
except:pass
# for i in [1,2,3,5,8,9,10,11,12]:
paths = ["05_08_A_01062017_Splice_shrink_FINISHED_LABEL_Bluecell_spinec91.ASC"]
syn_poses={}
syn_poses['05_08_A_01062017_Splice_shrink_FINISHED_LABEL_Bluecell_spinec91']=[(-5.56, -325.88, -451.42)]
for p in paths:
    print(p)
    cell_name=p[:-4]
    # dendogram = Dendogram('dend_only', p, add_sec2)#@#
    # dendogram.cumpute_distances(dendogram.cell.soma[0])#@#
    dendogram = Dendogram('dend_only', p, add_sec2,dots_loc=syn_poses[cell_name])
    dendogram.cumpute_distances(dendogram.cell.soma[0])
    max_y=dendogram.plot(save_folder2,title=save_folder2[:-2],ylabel="distance from soma (um)")

    dendogram = Dendogram('all', p, add_sec2, del_axon=False,dots_loc=syn_poses[cell_name])
    dendogram.cumpute_distances(dendogram.cell.soma[0])
    max_y = dendogram.plot(save_folder2,title=save_folder2[:-2],ylabel="distance from soma (um)")

    dendogram = Dendogram('dend_only_with_syn', p, add_sec,dots_loc=syn_poses[cell_name])
    dendogram.cumpute_distances(dendogram.cell.soma[0])
    max_y = dendogram.plot(save_folder,title=save_folder[:-2],ylabel="distance from soma (lamda)")

    dendogram = Dendogram('all_with_syn', p, add_sec, del_axon=False,dots_loc=syn_poses[cell_name])
    dendogram.cumpute_distances(dendogram.cell.soma[0])
    max_y = dendogram.plot(save_folder,title=save_folder[:-2],ylabel="distance from soma (lamda)")


