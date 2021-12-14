import pickle
import numpy as np
from glob import glob
#from correct_noise import clear_noise,fast_clear_noise
from tqdm import tqdm


def read_from_pickle(path, hz=False ,rest=False ):
	#print('opening '+path)
	with open(path, 'rb') as file:
		try:
			while True:
				object_file=pickle.load(file)
		except EOFError:
			pass
		file.close()
	return object_file
		# if hz and rest:
		# 	HZ = object_np[2]
		# 	Rest = object_np[1]
		# 	return data,HZ, Rest
		# elif rest:
		# 	Rest = object_np[1]
		# 	return data,Rest
		# elif hz:
		# 	HZ = object_np[1]
		# 	return data,HZ
		# else:
		#
		#

if __name__=='__main__':
	from correct_noise import clear_noise
	#pickle_path = '2017_05_08_A_4-5_stable_conc_aligned_selected_Moria.abf.pickle'
	#bring_phenomena(pickle_path)
	#pickle_path = '/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis/2017_05_08_A_0006.pickle'
	a,hz=read_from_pickle('/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/data/short_pulse/_8.p')
	# correct=read_phenomena('syn','/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/data',1)
	a=1
