import numpy as np
import pickle
import gc
from open_pickle import read_from_pickle
import matplotlib.pyplot as plt
from add_figure import add_figure
import signal
def SIGSEGV_signal_arises(signalNum, stack):
    print(f"{signalNum} : SIGSEGV arises")
signal.signal(signal.SIGSEGV, SIGSEGV_signal_arises)

def reshape_data(data):
	minlen = min(len(r) for r in data)
	new = []
	for i in data:
		if minlen - len(i) < 0:
			new.append(np.delete(i, np.s_[minlen - len(i):]))
		else:
			new.append(i)
	return new
def clear_phenomena_partial(phenomena,phenomena_name,part_name,base,std_max=3.3,start=None,end=None,correct_rest=False):
	temp=np.mean(phenomena,axis=0)
	add_figure('clear part of the graph by max','point','mV')
	for v in phenomena:
		plt.plot(v)
	# plt.plot(temp[start:end],'black',linewidth=7)
	plt.plot(np.arange(start, end), temp[start:end],'black',linewidth=7)
	plt.savefig(base+phenomena_name+'/place2clear_bymax_'+part_name)
	index2del=[]
	add_figure('clear noises from ' + phenomena_name, 'Vec_index', 'mV')
	count = 0
	for i, v in enumerate(phenomena):
		data_np=np.array(v[start:end])
		m = data_np.mean(axis=0)
		s = data_np.std(axis=0)
		plt.plot(v,alpha=0.2)
		if  max(data_np) > m+s*std_max or min(data_np) < np.mean(m)-s*std_max:
			plt.plot(v,'black',linewidth=1,label='clear by max')
			index2del.append(i)
			count+=1
		if correct_rest:
			phenomena[i]=v-m

	print(count,' is remove from '+part_name)
	plt.suptitle( str(count)+ ' remove from graph bymax')
	plt.savefig(base+phenomena_name+'/noise2clear_'+part_name+'_'+phenomena_name)
	plt.close()
	if correct_rest: return index2del,phenomena
	else:return index2del

def clear_phenomena(phenomena,phenomena_name,base,std_mean=3,std_max=6,bymax=False, bymean=True):
	index2del= []
	mean_V = np.mean(np.array(phenomena), axis=0)
	count1, count2 = 0, 0
	add_figure('clear noises from '+phenomena_name,'Vec_index','mV')
	for i, v in enumerate(phenomena):
		data_np=np.array(v)
		m = data_np.mean(axis=0)
		s = data_np.std(axis=0)
		plt.plot(v,alpha=0.2)
		if bymean:
			if np.mean(mean_V)>m+std_mean*s or np.mean(mean_V)<m-s*std_mean:
				plt.plot(v,'blue',alpha=0.3,linewidth=1,label='clear by mean')
				index2del.append(i)
				count1+=1

		if bymax:
			if  max(v) > np.mean(m)+s*std_max or min(v) < np.mean(m)-s*std_max:
				#phenomena.pop(i)
				plt.plot(v,'black',linewidth=1,label='clear by max')
				index2del.append(i)
				count2+=1
	if bymean and bymax:
		plt.suptitle(str(count1)+' remove by std_mean (blue) and '+str(count2)+ ' by std_max (black)')
	elif bymean and not bymax:
		plt.suptitle(str(count1)+' remove by std_mean - blue')
	elif bymax and not bymean:
		plt.suptitle( str(count2)+ ' remove by std_max - black')

	plt.savefig(base+phenomena_name+'/noise2clear_'+phenomena_name)
	plt.close()
	return index2del

def phenomena(t,T,base,x_units='S',Y_units='mV'):

	V,short_pulse,spike,syn,noise1,noise2,rest1,rest2,mean_V =[], [],[], [],[],[],[],[],[]
	for v in np.array(t):
		V.append(v)
		short_pulse.append(v[:10000])
		noise1.append(v[10000:22000])
		spike.append(v[22000:34000])
		syn.append(v[32000:40000])
		noise2.append(v[40000:])
		rest1.append(np.mean(v[10000:22000]))
		rest2.append(np.mean(v[40000:]))
		mean_V.append(np.mean(v))
	T_spike=T[0][22000:34000]
	T_syn=T[0][32000:40000]
	T_short_pulse=T[0][:10000]
	with open(base + '/V1/V.p', 'wb') as f:
		pickle.dump( [V,T], f)
	with open(base+'/syn/syn.p', 'wb') as f:
		pickle.dump( [np.array(syn)*t.units,T_syn], f)
	with open(base+'/short_pulse/short_pulse.p', 'wb') as f:
		pickle.dump([np.array(short_pulse)*t.units,T_short_pulse], f)
	with open(base + '/spike/spike.p', 'wb') as f:
		pickle.dump([np.array(spike)*t.units,T_spike], f)
	with open(base + '/noise1/noise1.p', 'wb') as f:
		pickle.dump(np.array(noise1), f)
	with open(base + '/noise2/noise2.p', 'wb') as f:
		pickle.dump(np.array(noise2), f)
	rest4list=np.mean([rest1,rest2],axis=0)
	REST=np.mean(rest4list)
	noise11 = [v - rest4list[i] for i, v in enumerate(noise1)]
	noise22 = [v - rest4list[i] for i, v in enumerate(noise2)]

	index2del_1=clear_phenomena(noise11,'noise1',base,bymax=True)
	index2del_2=clear_phenomena(noise22,'noise2',base,bymax=True)
	new_V=np.delete(V,list((index2del_1+index2del_2)) ,axis=0)- REST
	new_short_pulse=np.delete(short_pulse,list(set(index2del_1+index2del_2)) ,axis=0)
	new_spike=np.delete(spike,list(set(index2del_1+index2del_2)) ,axis=0)
	new_syn=np.delete(syn,list(set(index2del_1+index2del_2)),axis=0)
	new_rest4list=np.delete(rest4list,list(set(index2del_1+index2del_2)),axis=0)
	# duble_rest_short_pulse=np.mean(new_V[40000:])

	new_short_pulse = [v - new_rest4list[i]-np.mean((v-new_rest4list[i])[:3000]) for i, v in enumerate(new_short_pulse)]
	new_spike = [v - new_rest4list[i]-np.mean((v-new_rest4list[i])[6000:12000]) for i, v in enumerate(new_spike)]
	new_syn = [v - new_rest4list[i]-np.mean((v-new_rest4list[i])[:2500]) for i, v in enumerate(new_syn)]

	index2del_short_pulse1 = clear_phenomena(new_short_pulse,'short_pulse',base,std_mean=0.2)
	index2del_spike1 = clear_phenomena(new_spike,'spike',base,std_mean=0.2)
	index2del_syn1 = clear_phenomena(new_syn,'syn',base,std_mean=0.7)

	new_short_pulse1=np.delete(new_short_pulse,index2del_short_pulse1,axis=0 )
	new_syn1=np.delete(new_syn,index2del_syn1 ,axis=0)

	temp_short_pulse=np.mean(new_short_pulse1,axis=0)
	temp_syn=np.mean(new_syn1,axis=0)
	short_pulse_time2clear1=np.argmax(temp_short_pulse)-10
	short_pulse_time2clear2=np.argmin(temp_short_pulse)+700
	syn_time2clear1=np.argmax(temp_syn)-100
	syn_time2clear2=np.argmax(temp_syn)+200
#@# add a pickle t save this places
	index2del_short_pulse_begining,new_short_pulse1 = clear_phenomena_partial(new_short_pulse1, 'short_pulse','begining', base ,start=short_pulse_time2clear1-500,end=short_pulse_time2clear1,correct_rest=True)
	index2del_short_pulse_end = clear_phenomena_partial(new_short_pulse1, 'short_pulse','end', base ,start=short_pulse_time2clear2,end=short_pulse_time2clear2+1000)
	index2del_syn_begining,new_syn1 = clear_phenomena_partial(new_syn1, 'syn','begining', base, start=syn_time2clear1-500,end=syn_time2clear1,correct_rest=True)
	new_syn1 = np.delete(new_syn1, list(set(index2del_syn_begining)), axis=0)
	new_syn11=[]
	for v in new_syn1:
		v = v - np.mean(v[syn_time2clear1 - 500:syn_time2clear1])
		new_syn11.append(v)
	plt.plot(v[syn_time2clear1 - 500:syn_time2clear2 + 1000])

	index2del_syn_begining2,new_syn2 = clear_phenomena_partial(new_syn1, 'syn','begining2', base, start=syn_time2clear1-700,end=syn_time2clear1,correct_rest=True)
	index2del_syn_end = clear_phenomena_partial(new_syn2, 'syn','end', base, start=syn_time2clear2,end=syn_time2clear2+1000,std_max=3.3)
	index2del_syn_end2 = clear_phenomena_partial(new_syn2, 'syn','end2', base, start=syn_time2clear2+1000,end=syn_time2clear2+3000,std_max=4)


	new_short_pulse2 = np.delete(new_short_pulse1, list(set(index2del_short_pulse_begining+index2del_short_pulse_end)), axis=0)+ REST
	new_spike2 = np.delete(new_spike, index2del_spike1, axis=0) + REST
	new_syn2 = np.delete(new_syn1, list(set(index2del_syn_begining+index2del_syn_begining2+index2del_syn_end+index2del_syn_end2)), axis=0)+ REST

	names=['short_pulse','spike','syn']
	for i,phenomena in enumerate([new_short_pulse2,new_spike2,new_syn2]):
		plt.close()
		add_figure('clear '+names[i],'index',t.units)
		for v in phenomena:
			plt.plot(v,alpha=0.1,lw=0.5,color='grey')
		plt.plot(np.mean(phenomena,axis=0),'black',lw=3)
		plt.savefig(base+names[i]+'/clear_'+names[i])
		with open(base +names[i]+'/clear_'+names[i]+'.p', 'wb') as f:
			pickle.dump( [np.array(phenomena),eval('T_'+names[i])], f)

	with open(base + '/V1/clear_V.p', 'wb') as f:
		pickle.dump( np.array(new_V), f)

	for i,phenomena in enumerate([new_short_pulse2,new_spike2,new_syn2]):
		add_figure('mean '+names[i],eval('T_'+names[i])[0].units,t.units)
		mean=np.mean(phenomena,axis=0)
		plt.plot(eval('T_'+names[i]),mean)
		plt.savefig(base+names[i]+'/mean_'+names[i])
		with open(base +names[i]+'/mean_'+names[i]+'.p', 'wb') as f:
			pickle.dump( [mean*t.units,eval('T_'+names[i])], f)
	with open(base + '/V1/mean_V.p', 'wb') as f:
		pickle.dump( [np.mean(np.array(new_V),axis=0)*t.units,T[0]], f)

	add_figure('mean short syn',T_syn.units,t.units)
	plt.plot(T_syn[syn_time2clear1-500:syn_time2clear2+1000],np.mean(new_syn2,axis=0)[syn_time2clear1-500:syn_time2clear2+1000])
	plt.close()
	with open(base +'syn'+'/mean_syn_short.p', 'wb') as f:
		pickle.dump([np.mean(new_syn2,axis=0)[syn_time2clear1-500:syn_time2clear2+1000] * t.units,T_syn[syn_time2clear1-500:syn_time2clear2+1000]], f)

	E_pas_syn = np.mean([new_short_pulse2[i][short_pulse_time2clear1 - 500:short_pulse_time2clear1] for i in range(len(new_short_pulse2))])
	E_pas_short_pulse=np.mean([new_syn2[i][syn_time2clear1-500:syn_time2clear1] for i in range(len(new_syn2))])
	E_pases=[E_pas_short_pulse,E_pas_syn]
	point2calculate_E_pas=[short_pulse_time2clear1,syn_time2clear1]
	names2=['short_pulse','syn']
	for i, phenomena in enumerate([new_short_pulse2, new_syn2]):
		mean=np.mean(phenomena,axis=0)
		with open(base +names2[i]+'/mean_'+names2[i]+'_with_parameters.p', 'wb') as f:
			pickle.dump({'mean':[mean * t.units, eval('T_' + names2[i])],'E_pas':E_pases[i],'points2calsulate_E_pas':point2calculate_E_pas }, f)
	#add to the other currents for I-V curve
	add_figure('I_V curve_together', 'points', t.units)
	plt.plot(new_short_pulse2)
	'/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/data/traces_img/2017_05_08_A_0006/-50pA.png'
	plt.savefig('data/traces_img/2017_05_08_A_0006/-50pA')
	with open('data/traces_img/2017_05_08_A_0006/-50pA.p', 'wb') as f:
		pickle.dump({'mean': [np.mean(new_short_pulse2,axis=0) * t.units, T_short_pulse], 'E_pas': E_pases[i],}, f)
	a=1
	return REST,np.mean(new_short_pulse2,axis=0)* t.units,T_short_pulse


def one_data( t,T, base, cut_part=25000):
	idx = np.where(t > 10)[0][0] + cut_part
	t1 = t[idx:]
	t2=t1
	V = []
	short_pulse = []
	spike = []
	syn = []
	noise1=[]
	noise2 = []

	while True: #split 1 picture to spyke,synapse,short_pulse and noise
		try:
			idx = np.where(t2 > 10)[0]
			if len(idx)==0:
				break
			idx=idx[0]+cut_part
			V.append(t2[:idx])
			short_pulse.append(t2[11000:17000])
			syn.append(t2[34000:36500])
			spike.append(t2[44000:47000])
			#noise1.append(t2[17000:34000])
			noise2.append(t2[47000:])
			t2 = t2[idx:]

		except:
			break

		syn = reshape_data(syn)
#	V = reshape_data(V)
	with open(base+'V.p', 'wb') as f:
		pickle.dump( np.array(V),f)
	with open(base+'syn.p', 'wb') as f:
		pickle.dump( np.array(syn),f)
	with open(base+'short_pulse.p', 'wb') as f:
		pickle.dump(np.array(short_pulse), f)
	with open(base + 'spike.p', 'wb') as f:
		pickle.dump( np.array(spike), f)
	with open(base + 'noise.p', 'wb') as f:
		pickle.dump( np.array(noise2), f)
	del t2,t1,t,syn,spike,noise2,short_pulse
	gc.collect()


if __name__=='__main__':
	path='/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/data/traces_img/2017_05_08_A_4-5_stable_conc_aligned_selected_Moria/first_channel.p'
	t,T=read_from_pickle(path)
	phenomena(t,T, 'data/')
	b=1