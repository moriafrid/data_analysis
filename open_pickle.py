import pickle
import numpy as np
from matplotlib import pyplot as plt
from glob import glob
from correct_noise import clear_noise,fast_clear_noise

def read_from_pickle(path, hz=False ,rest=False ):
	#print('opening '+path)
	with open(path, 'rb') as file:
		try:
			while True:
				object_file=pickle.load(file)
		except EOFError:
			pass
		object_np=np.squeeze(np.array(tuple(object_file.values())))
		file.close()
		return object_np

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
	for path in glob(base+phenomena+'/*.p'):
		data.append(read_from_pickle(path))
		#data.append(pickle.load(open(path,'rb')))
	correct=fast_clear_noise(data, str(phenomena), base, jumps)
#	correct=clear_noise(data, str(phenomena), base, jumps)
	return correct

if __name__=='__main__':
	from correct_noise import clear_noise
	#pickle_path = '2017_05_08_A_4-5_stable_conc_aligned_selected_Moria.abf.pickle'
	#bring_phenomena(pickle_path)
	#pickle_path = '/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis/2017_05_08_A_0006.pickle'
	correct=read_phenomena('syn','/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis/data/',100)
