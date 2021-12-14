import numpy as np
import pickle
import gc

def reshape_data(data):
	minlen = min(len(r) for r in data)
	new = []
	for i in data:
		if minlen - len(i) < 0:
			new.append(np.delete(i, np.s_[minlen - len(i):]))
		else:
			new.append(i)
	return new
def split_with_generator(t2,size):
	yield t2[size]

def one_data(index, t, hz, base, cut_part=25000):
	idx = np.where(t > 10)[0][0] + cut_part
	t1 = t[idx:]
	t2=t1
	V = []
	short_pulse = []
	spike = []
	syn = []
	noise1 = []
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
	with open(base+'syn/_'+str(index)+'.p', 'wb') as f:
		pickle.dump({'syn_'+str(index): np.array(syn),'HZ':hz},f)
	with open(base+'short_pulse/_'+str(index)+'.p', 'wb') as f:
		pickle.dump({'short_pulse_'+str(index): np.array(short_pulse),'HZ':hz}, f)
	with open(base + 'spike/_' + str(index) + '.p', 'wb') as f:
		pickle.dump({'spike_' + str(index): np.array(spike),'HZ':hz}, f)
	with open(base + 'noise/_' + str(index) + '.p', 'wb') as f:
		pickle.dump({'noise2_' + str(index): np.array(noise2),'HZ':hz}, f)
	del t2,t1,t,syn,spike,noise2,short_pulse,idx
	gc.collect()


	#	'noise1_'+str(index+i): np.array(noise1)),
	#}, open(base+'noise/1_'+str(index+i)+'.p', 'wb'))
	#pickle.dump({


	#return V,short_pulse,syn,spike,noise1,noise2
	#return short_pulse,syn,spike,noise2

