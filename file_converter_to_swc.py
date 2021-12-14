from neuron import h, gui
import numpy as np

folder_='/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/'
h.load_file("import3d.hoc")
h.load_file("nrngui.hoc")
h.load_file('stdlib.hoc')
h.load_file("stdgui.hoc")

def run(id, prev_id,sec,type, print_=True):
    sec_points = np.array([list(i) for i in sec.psection()['morphology']['pts3d']])
    if print_:
        print(sec.name(), len(sec.children()))
    for point in sec_points:
        swc_file.write(str(id)+' '+str(type)+' '+
                       ' '.join(point[:3].round(4).astype(str).tolist()) +
                       ' ' + str(round(point[3]/ 2.0, 4))+' '+str(prev_id)+'\n')
        prev_id=id
        id+=1
    for child in sec.children():
        id=run(id,prev_id,child, type, print_=print_)
    return id

class Cell: pass
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

fname = folder_+"05_08_A_01062017_Splice_shrink_FINISHED_LABEL_Bluecell_spinec91.ASC"
cell=mkcell(fname)
sp2 = h.PlotShape()
sp2.color_all(3)
sp2.show(0)  # show diameters
soma_points = np.array([list(i) for i in cell.soma[0].psection()['morphology']['pts3d']]).mean(axis=0)
swc_file = open(folder_+'try1.swc', 'w')
swc_file.write('# generated by Vaa3D Plugin sort_neuron_swc\n')
swc_file.write('# source file(s): '+fname+'\n')
swc_file.write('# id,type,x,y,z,r,pid\n')
swc_file.write('1 1 '+
               ' '.join(soma_points[:3].round(4).astype(str).tolist())+
               ' '+str(round(cell.soma[0].diam/2.0, 4))+' -1\n')

id=2
for child in cell.soma[0].children():
    if child in cell.dend:
        type=3 #2 for dend
    elif child in cell.axon:
        type=2
    else:
        raise Exception('no type chosen')
    id=run(id,1,child,type, print_=type==2)

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

cell=instantiate_swc('/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/try1.swc')
print (cell)
sp = h.PlotShape()
sp.show(0)  # show diameters
a=1