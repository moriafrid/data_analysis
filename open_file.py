#! /usr/bin/env python3

#import sys
#sys.path.append("/home/hines/lib64/python")
#import neuron
#from neuron import h

import os
from neo import io
from glob import glob
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import numpy as np
import  os
import pickle


abf_files=[]
#morph_file_path = 
#os.path.join("C:\\Users\\moria.fridman\\Documents\\Segev_Lab\\project\\nueron1\\dropbox", 
#"morphology_correctedZ.hoc")
#h.load_file(morph_file_path)

folder_="."
base = "data2"
try: os.mkdir('data2')
except:pass
try:os.mkdir(base + '/traces_img')
except:pass
#os.chdir("dropbox/2017_05_08_A_0006.abf")
print('start2')
for f in os.listdir(folder_):
    if f.endswith(".abf"):
        print(f)
        dicty = {}
        abf_files.append(f)
        save_folder = os.path.join(base,'traces_img', f[f.rfind('/')+1:-4])
        try:
            #os.mkdir(folder_+'abf_naw project folder inside the auto-syncedme.abf')
            os.mkdir(save_folder)
        except Exception as ex:
            print(str(ex.args[0]) + " save folder = " + save_folder)
            pass
        r=io.AxonIO(f)
        bl = r.read_block(lazy=False)
        for i in range(len(bl.segments)):
            print(i)
            segment = bl.segments[i]
            hz = segment.analogsignals[0].sampling_rate
            t = [np.array(segment.analogsignals[0]) for segment in bl.segments]
            t = np.array(t).flatten()
            T = np.arange(0, len(t), 1) * (1000.0 / hz)
            dicty[str(f)+'/_'+str(i)]=[t,T]

            plt.close()
            plt.plot(T, t)
            plt.savefig(save_folder+'/_'+str(i))
            plt.close()
        with open(str(f) + '.pickle', 'wb') as g:
            pickle.dump(dicty, g, protocol=pickle.HIGHEST_PROTOCOL)
dicty2={'abf_files':abf_files}
with open('abf_name.pickle', 'wb') as g:
    pickle.dump(dicty2, g)
