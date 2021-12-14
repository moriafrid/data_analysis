import numpy as np
from neuron import h,gui
import signal
from glob import glob

def synaptic_loc(cell,syn_poses,del_axon=True):
    print('unsure axon is deleted before synaptic_loc.py or del_axon=False')
    #syn_poses should be (x,y,z) coordinates
    h.load_file("import3d.hoc")
    h.load_file("nrngui.hoc")
    # file_name=glob(cell_name+'*.ASC')
    # cell = mkcell(file_name[0])
    secs,dends,dists,dends_name=[],[],[],[]
    for i in range(len(syn_poses)):
        secs.append(None)
        dends.append(None)
        dists.append(10000)
        dends_name.append(None)
    for sec in h.allsec():
        if del_axon==False:
            if sec in cell.axon: continue
        lens = []
        initial_point = np.array([sec.x3d(0), sec.y3d(0), sec.z3d(0)])
        for i in range(sec.n3d()):
            lens.append(np.linalg.norm(initial_point - np.array([sec.x3d(i), sec.y3d(i), sec.z3d(i)])))
            initial_point = np.array([sec.x3d(i), sec.y3d(i), sec.z3d(i)])
        total_len = np.sum(lens)
        accumalate_len = 0
        initial_point = np.array([sec.x3d(0), sec.y3d(0), sec.z3d(0)])
        for i in range(sec.n3d()):
            dend_pos = np.array([sec.x3d(i), sec.y3d(i), sec.z3d(i)])
            accumalate_len += np.linalg.norm(initial_point - dend_pos)
            initial_point = dend_pos
            for j, syn_pos in enumerate(syn_poses):
                if np.linalg.norm(syn_pos - dend_pos) < dists[j]:
                    dists[j] = np.linalg.norm(syn_pos - dend_pos)
                    secs[j] = [sec, accumalate_len / total_len]
                    dends[j]=[sec,round(accumalate_len / total_len,3)]
                    dends_name[j]=[str(sec)[str(sec).find('>')+2:],round(accumalate_len / total_len,3)]
    return {'place_name':dends_name,'place_as_sec':dends}
def syn_dis_from_soma(cell,synapses_locations):
    h.distance(0, 0.5, sec=cell.soma[0])
    synapses_dis_from_soma = []
    for syn_loc in synapses_locations:
        synapses_dis_from_soma.append([syn_loc[0], h.distance(eval('cell.' + syn_loc[0])(syn_loc[1]))])
    return synapses_dis_from_soma
def SIGSEGV_signal_arises(signalNum, stack):
    print(f"{signalNum} : SIGSEGV arises")

class Cell(): pass
def mkcell(fname):
    # def to read ACS file
    loader = h.Import3d_GUI(None)
    loader.box.unmap()
    loader.readfile(fname)
    c = Cell()
    loader.instantiate(c)
    return c
if __name__=='__main__':


    signal.signal(signal.SIGSEGV, SIGSEGV_signal_arises)

    cell_name = "170830HuSHS2C1IN0toIN3"
    #cell = mkcell(fname)
    syn_poses={}
    syn_poses["170830HuSHS2C1IN0toIN3"] = [np.array([198.04, 51.55, 11.31]),
                 np.array([87.81, -41.95, -13.15]),
                 np.array([-143.67, 210.14, 23.0])]
    cell_name='05_08_A_01062017_Splice_shrink_FINISHED_LABEL_Bluecell_spinec91'
    syn_poses['05_08_A_01062017_Splice_shrink_FINISHED_LABEL_Bluecell_spinec91'] = [(-5.56, -325.88, -451.42)]
    synapses_locations=synaptic_loc(cell_name,syn_poses[cell_name])
    a=1