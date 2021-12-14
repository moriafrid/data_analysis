import pickle
import numpy as np
from glob import glob
from correct_noise import clear_noise,fast_clear_noise
from tqdm import tqdm


def read_from_pickle(path, hz=False ,rest=False ):
	#print('opening '+path)
	with open(path, 'rb') as file:
		try:
			while True:
				object_file=pickle.load(file)
		except EOFError:
			pass
		return object_file


def bring_phenomena(pickle_path):
	all_split=read_from_pickle(pickle_path)
	short_pulse=all_split['short_pulse']
	spike=all_split['spike']
	noise1=all_split['noise1']
	noise2=all_split['noise2']
	syn=all_split['syn']
	V=all_split['V']
	return V,spike,syn,short_pulse,noise1,noise2

def read_phenomena(phenomena,base,jumps):
	data=[]
	hz=[]
	x_data=[]
	for path in tqdm(glob(base+phenomena+'/*.p')):
		open_data,open_hz=read_from_pickle(path, hz=True)
		data.append(open_data)
		x_data.append(np.arange(0, 9992) + np.zeros((len(open_data),1)))
		#data.append(pickle.load(open(path,'rb')))

	correct,mean=fast_clear_noise(data,hz, str(phenomena), base, jumps)
#	correct=clear_noise(data, str(phenomena), base, jumps)
	return correct,mean,hz

if __name__=='__main__':
	from correct_noise import clear_noise
	#pickle_path = '2017_05_08_A_4-5_stable_conc_aligned_selected_Moria.abf.pickle'
	#bring_phenomena(pickle_path)
	#pickle_path = '/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis/2017_05_08_A_0006.pickle'
	a,hz=read_from_pickle('/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/data/short_pulse/_8.p')
	correct=read_phenomena('syn','/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/data',1)
	a=1
