import os
from neo import io
from glob import glob
from matplotlib import pyplot as plt
import numpy as np
from tqdm import tqdm
import pickle
from open_one_data import one_data

def split2phenomena(folder_):
	try:	os.mkdir(folder_ + 'data')
	except:	pass
	base = folder_+'data/'
	for phen in ['V1', 'short_pulse', 'syn', 'spike', 'noise']:
		try:os.mkdir(base + phen)
		except:pass
	abf_files = glob(folder_ + '*.abf')


	for f in abf_files[:]: #take the abf file (from 3)
		print(f)
		save_folder = base+'traces_img/'+f[f.rfind('/')+1:-4]
		try: os.mkdir(save_folder)
		except:pass

		r=io.AxonIO(f)
		bl = r.read_block(lazy=False)
		for i in tqdm(range(len(bl.segments))): # on all the picture (from 1+92+10)
			segment =bl.segments[i]
			hz = segment.analogsignals[0].sampling_rate
			t = [np.array(segment.analogsignals[0]) for segment in bl.segments]
			t = np.array(t).flatten()
			T = np.arange(0, len(t), 1)*(1000.0/hz)
			plt.close()
			plt.plot(T, t)
			plt.savefig(save_folder+'/'+'_'+str(i))
			plt.close('all')

			#split to syn, short_pulse, spike ,noise
			if f=='/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis/2017_05_08_A_4-5_stable_conc_aligned_selected_Moria.abf':
				one_data(i,t,hz,base)
			else:
				with open(save_folder+'/_' + str(i) + '.p', 'wb') as f:
					pickle.dump({str(i): [np.array(t),np.array(T)]},f)




if __name__=='__main__':
	folder_ = '/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis/'
	split2phenomena(folder_)




