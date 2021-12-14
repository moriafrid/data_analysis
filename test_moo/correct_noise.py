import sys,os
file_dir = os.path.dirname(__file__)
sys.path.append(file_dir)
import matplotlib
from neo import io
from matplotlib import pyplot as plt
import numpy as np
from open_one_data import one_data
from add_figure import add_figure
matplotlib.use('agg')
import pickle
from tqdm import tqdm
#hendel function
def reshape_data(data):
	minlen = min(len(r) for r in data)
	new = []
	for i in data:
		if minlen - len(i) < 0:
			new.append(np.delete(i, np.s_[minlen - len(i):]))
		else:
			new.append(i)
	return new

##long function that show the levels of clearing the noise:
def clear_noise(data,name,base,jumps):
	#save in folder

	save_folder = os.path.join(base,name)
	try:
		os.mkdir(save_folder)
	except FileExistsError as ex:
		print(str(ex.args[0]) + " save f￼older = " + save_folder)
	#the initial data before noise clearning
	add_figure(name,'time [ms]','V [mV]')
	for v in data:
		rest = np.mean(v[:10000])
		plt.plot(v[::jumps]-rest, alpha=0.5)
	plt.savefig(save_folder+'/befor_noise_remove')
	plt.close()
	print('the initual data is printing')

	data2=[]
	data=reshape_data(data)
	#
	for i in data:
		data2.append(i-np.mean(i[:500]))
	data_np = np.array(data2)
	## plot the pulse with std and the mean
	plt.close()
	plt.figure()
	for i in data2:
		plt.plot(i[::jumps],color='gray', alpha = 0.05)
	m=data_np.mean(axis=0)[::jumps]
	s=data_np.std(axis=0)[::jumps]
	plt.fill_between(  np.arange(len(m+s)),m+s, color='green', alpha=0.5)
	plt.fill_between( np.arange(len(m-s)), m-s, color='green', alpha=0.5)
	plt.fill_between( np.arange(len(m+s)), m+s*2, color='green', alpha=0.1)
	plt.fill_between( np.arange(len(m-s)), m-s*2, color='green', alpha=0.1)
	plt.savefig(save_folder+'/std_fill')
	plt.close()
	print('the std is printing')
	plt.figure()
	for i in data_np:
		plt.plot(i[::jumps] - m, color='k', alpha=0.2)
		plt.plot(+s, color='r')
		plt.plot(s * 2, color='y')
		plt.plot(-s, color='r')
		plt.plot(-s * 2, color='y')
	plt.savefig(save_folder + '/std_flat_data')
	plt.close()

	# remove the mosts noise pulses
	correct = []
	correct_num = 0
	wrong_num=0
	plt.figure()
	for i in data_np:
		b = i[::10] - m
		if b[np.argmax(b)] < 2 * s[np.argmax(b)] and b[np.argmin(b)] > -2 * s[np.argmin(b)]:
			plt.plot(b, color='k', alpha=0.2)
			correct.append(i)
			correct_num += 1
		else:
			wrong_num+=1
		plt.plot(+s, color='r')
		plt.plot(s * 2, color='y')
		plt.plot(-s, color='r')
		plt.plot(-s * 2, color='y')
	# plot th# plot the pulse with std and the mean
	plt.savefig(save_folder + '/after_clean_noise_flat')
	plt.close()
	print('in '+name+' the number of right traces is '+str(correct_num)+'/'+str(correct_num+wrong_num))


	plt.figure()
	for i in correct:
		plt.plot(i[::jumps], alpha=0.2)
	plt.savefig(save_folder + '/data_after_clean_noise')
	plt.close()

##short function to clear the noises
def fast_clear_noise(data,HZ,name,base,jumps):
	print('start clear noise from'+name)
	save_folder = os.path.join(base,name)
	#the initial data before noise clearning
	data2=[]
	data=reshape_data(data)
	data = np.array(data)
	a,b,c=np.shape(data)
	data=data.reshape(a*b,c)
	add_figure(name+' before clean noise', 'time [ms]', 'V [mV]')
	rest_list=[]
	for v in data:
		rest = np.mean(v[:10000])
		plt.plot(np.arange(len(v[::jumps]))*(1000.0/hz),v[::jumps]-rest, alpha=0.5)
		data2.append(v - np.mean(v[:500]))
		rest_list.append(rest)
	plt.savefig(save_folder+'/befor_noise_remove')
	plt.close()

	data_np=np.array(data2)
	# remove the mosts noise pulses and count them
	correct = []
	correct_num = 0
	wrong_num=0
	add_figure(name,'time [ms]','V [mV]')
	#a,b,c=np.shape(data_np)
	#data_np=data_np.reshape(a*b,c)
	m=data_np.mean(axis=0)[::jumps]
	s=data_np.std(axis=0)[::jumps]
	for i in data_np:
		b = i[::jumps] - m
		if b[np.argmax(b)] < 2 * s[np.argmax(b)] and b[np.argmin(b)] > -2 * s[np.argmin(b)]:
			correct.append(i)
			correct_num += 1
		else:
			wrong_num+=1

	add_figure(name+' after clean noise', 'time [ms]', 'V [mV]')
	for i in tqdm(enumerate(correct,HZ)):
		plt.plot(range(len(i[::jumps]))*(1000.0/hz),i[::jumps], alpha=0.2)
	mean=np.mean(correct,axis=0)
	plt.plot(range(len(mean[::jumps])),mean[::jumps],color='black')
	plt.savefig(save_folder + '/after_clean_noise')
	plt.close()
	print('in '+name+' the number of right traces is '+str(correct_num)+'/'+str(correct_num+wrong_num))
	add_figure(name+' mean', 'time [ms]', 'V [mV]')
	plt.plot(range(len(i[::jumps]))*(1000.0/hz),mean[::jumps])
	plt.savefig(save_folder + '/mean')
	plt.close()
	#save the corrects in pickle
	with open('data/'+name+'_correct.p', 'wb') as f:
		pickle.dump({		name+'_correct': correct, 'rest':rest_list, 'HZ': hz}, f)
	with open('data/'+name+'_mean.p', 'wb') as f:
		pickle.dump({		name+'_correct': correct, 'rest':mean(rest_list), 'HZ': hz}, f)

	return correct, mean



if __name__=='__main__':
	# one data
	base='/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis/data2'
	f='/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis/2017_05_08_A_4-5_stable_conc_aligned_selected_Moria.abf'
	V, short_pulse, syn, spike, noise1, noise2 = one_data(5, f)
	clear_noise(V, 'V1', base, 10)
	clear_noise(syn, 'syn1', base, 10)
	clear_noise(short_pulse, 'short_pulse1', base, 10)
	clear_noise(spike, 'spike1', base, 10)

	pickle_path='/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis/all_split.p'
	clear_noise(V, 'V', base, 20)
	clear_noise(syn, 'syn', base, 20)
	clear_noise(short_pulse, 'short_pulse', base, 20)
	clear_noise(spike, 'spike', base, 20)

