from neo import io
from glob import glob
from matplotlib import pyplot as plt
import numpy as np
from tqdm import tqdm
import pickle
from open_one_data import phenomena
import os
from add_figure import add_figure
import quantities as pq
from IV_curve import I_V_curve,sepereat_by_current,find_maxi
from check_dynamics import check_dynamics

def split2phenomena(folder_):
	try:	os.mkdir(folder_ + 'data')
	except FileExistsError:	pass
	base = folder_+'data/'
	try:os.mkdir(base+'traces_img')
	except FileExistsError:pass
	for phen in ['V1', 'short_pulse', 'syn', 'spike', 'noise']:
		try:os.mkdir(base + phen)
		except FileExistsError:pass
	abf_files = glob(folder_ + '*.abf')

	for f in abf_files[:]: #take the abf file (from 3)
		print(f,flush=True)
		save_folder = base+'traces_img/'+f[f.rfind('/')+1:-4]
		try: os.mkdir(save_folder)
		except FileExistsError:pass

		r=io.AxonIO(f)
		bl = r.read_block(lazy=False)
		hz = [np.array(segment.analogsignals[0].sampling_rate) for segment in bl.segments]
		t,T,t1,t2=[],[],[],[]
		second_channel=True
		for segment in tqdm(bl.segments):
			t_i = segment.analogsignals[0]
			channel1 = [v[0] for v in np.array(t_i)]
			t.append(t_i)
			t1.append(channel1)

			if second_channel:
				channel2 = [v[1] for v in np.array(t_i)]
				t2.append(channel2)
			T_i = np.linspace(segment.analogsignals[0].t_start, segment.analogsignals[0].t_stop, int(len(t_i)))
			T.append(T_i)

		with open(save_folder + '/full_channel.p', 'wb') as fr:
			pickle.dump( [ t,T], fr)
		with open(save_folder + '/first_channel.p', 'wb') as fr:
			pickle.dump( [ np.array(t1)*t_i.units,T], fr)
		add_figure(f[f.rfind('/') + 1:-4]+'\n first_channel',T[0].units,t[0].units)
		plt.plot(np.array(T).flatten(),np.array(t1).flatten())
		plt.savefig(save_folder+ '/first_channel.png')
		plt.savefig(save_folder+ '/first_channel.pdf')

		if second_channel:
			with open(save_folder + '/second_channel.p', 'wb') as fr:
				pickle.dump([ np.array(t2)*t_i.units,T], fr)
			plt.close()
			add_figure(f[f.rfind('/') + 1:-4]+'\n second_channel',T[0].units,t[0].units)
			plt.plot(np.array(T).flatten(), np.array(t2).flatten())
			plt.savefig(save_folder+ '/second_channel.png')
			plt.savefig(save_folder+ '/second_channel.pdf')

		#split to syn, short_pulse, spike ,noise
		if f==folder_+'2017_05_08_A_4-5_stable_conc_aligned_selected_Moria.abf':
			print(f,'correct one_data')
			REST,short_pulse,T_short_pulse=phenomena(np.array(t1)*t_i.units,T,base,x_units=T[0].units,Y_units=t_i.units)
		if f==folder_+'2017_05_08_A_0006.abf':
			f_name='2017_05_08_A_0006'
			save_folder_IV_curve=folder_+'data/traces_img/'+f_name+'/'
			I = [-200, -160, -120, -80, -40, -0, 40, 80, 120, 160]
			# print(f,'correct IV_curve')
			maxi=sepereat_by_current(np.array(t1)*t_i.units,T,I,f_name)
			maxi = np.append(maxi, find_maxi(np.array(short_pulse)-REST, f_name))
			I.append(-50)
			with open(save_folder_IV_curve + 'max_vol_curr_inj.p', 'wb') as fr:
				pickle.dump([maxi*short_pulse.units,I*pq.pA], fr)
			I_V_curve(maxi,I*pq.pA,save_folder_IV_curve)
			check_dynamics(short_pulse, T_short_pulse, folder_+'data/')

if __name__=='__main__':
    folder_ = '/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/'
    split2phenomena(folder_)




