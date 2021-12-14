from glob import glob
import sys,os
file_dir = os.path.dirname(__file__)
sys.path.append(file_dir)
from matplotlib import pyplot as plt
import numpy as np
from add_figure import add_figure
import pickle
from tqdm import tqdm
from open_pickle import read_from_pickle
import time

def reshape_data(data):
	minlen = min(len(r) for r in data)
	new = []
	for i in data:
		if minlen - len(i) < 0:
			new.append(np.delete(i, np.s_[minlen - len(i):]))
		else:
			new.append(i)
	return new

def clear_noise(phenomena,base,STDtime2clear=2,jumps=1):
	name=str(phenomena)
	print('start clear noise from '+name,flush=True)
	save_folder = os.path.join(base, name)

	data=[]
	x_data=[]
	hz_list=[]
	for path in tqdm(glob(base + phenomena + '/*.p')):
		open_data, open_hz = read_from_pickle(path, hz=True)
		data.append(open_data)
		hz=np.array(open_hz)
		x_data.append(np.arange(0, np.shape(open_data)[1])*1000/hz + np.zeros((len(open_data), 1)))
		hz_list.append(hz)
		#x_data.append(np.arange(0, 9992) + np.zeros((len(open_data),1)))
	#	the initial data before noise clearning

	data=reshape_data(data)
	x_data=reshape_data(x_data)
	data = np.array(data)
	x_data = np.array(x_data)
	a,b,c=np.shape(data)
	data=data.reshape(a*b,c)
	x_data=x_data.reshape(a*b,c)
	add_figure(name+' before clean noise', 'time [ms]', 'V [mV]')

	data2 = []
	rest_list=[]
	print('plot the traces of '+name+ ' before the noise removal',flush=True)
	for i,v in enumerate(tqdm(data)):
		rest = np.mean(v[:900])
		plt.plot(x_data[i][::jumps+3],v[::jumps+3]-rest, alpha=0.5)
		data2.append(v - np.mean(v[:900]))
		rest_list.append(rest)
	start_program = time.time()
	plt.savefig(save_folder+'/befor_noise_remove.png')
	#plt.savefig(save_folder+'/befor_noise_remove.pdf')

	end_program = time.time()
	print('the time taken for savefig before noise removal = ',end_program-start_program)
	plt.close()

	data_np=np.array(data2)
	# remove the mosts noise pulses and count them
	correct = []
	x_correct=[]
	correct_num = 0
	wrong_num=0
	add_figure(name,'time [ms]','V [mV]')
	m=data_np.mean(axis=0)
	s=data_np.std(axis=0)
	print('remove the wrong traces of '+name,flush=True)

	for i,trace in enumerate(tqdm(data_np)):
		b = trace - m
		if b[np.argmax(b)] < STDtime2clear * s[np.argmax(b)] and b[np.argmin(b)] > -STDtime2clear * s[np.argmin(b)]:
			correct.append(trace)
			x_correct.append(x_data[i])
			correct_num += 1
		else:
			wrong_num+=1

	add_figure(name+' after clean noise', 'time [ms]', 'V [mV]')
	print('plot the traces of '+name+ ' after the noise removal', flush= True)
	for i,v in enumerate(tqdm(correct)):
		plt.plot(x_correct[i][::jumps],v[::jumps], alpha=0.2)
	mean=np.mean(correct,axis=0)
	x_mean=np.mean(x_correct,axis=0)
	plt.plot(x_mean,mean,color='black')
	start_program=time.time()
	plt.savefig(save_folder + '/after_clean_noise')
	plt.savefig(save_folder+'/after_clean_noise.pdf')

	end_program=time.time()
	print('the time taken for savefig after noise removal = ', end_program - start_program)
	plt.close()

	print('in '+name+' the number of right traces is '+str(correct_num)+'/'+str(correct_num+wrong_num))
	add_figure(name+' mean', 'time [ms]', 'V [mV]')
	plt.plot(x_mean,mean,color='black')
	plt.savefig(save_folder + '/mean')
	plt.savefig(save_folder + '/mean.pdf')

	plt.close()

	#save the corrects in pickle
	with open('data/'+name+'_correct.p', 'wb') as f:
		pickle.dump({		name+'_correct': [correct,x_correct], 'rest':rest_list, 'HZ': np.mean(hz_list)}, f)
	with open('data/'+name+'_mean.p', 'wb') as f:
		pickle.dump({		name+'_mean': [mean,x_mean], 'rest':np.mean(rest_list), 'HZ': np.mean(hz_list)}, f)
	if np.std(rest_list)*80>abs(np.mean(rest_list)):
		print('the rest list standart deviation is too large and it cost ',np.std(rest_list),' from ',np.mean(rest_list),flush=True)
	if np.mean(hz_list)<np.std(hz_list)*1000:
		print("there is diffrence betwwen the measurments's hz", flash=True)
		print('the hz list is ',hz_list,flush=True)
	return mean ,x_mean ,correct ,x_correct

if __name__=='__main__':
	STDtime2clear=2
	folder_data = '/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/data/'
	correct, mean, x_correct, x_mean=	clear_noise('short_pulse', folder_data,STDtime2clear=2.5)
	a=1

