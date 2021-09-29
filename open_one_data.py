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
def one_data(index,t,hz,base):
	idx = np.where(t > -10)[0][0] + 25000
	t1 = t[idx:]
	t1 = t1.reshape(2, len(t1) // 2)
	for i in (range(1)):
		t2=t1[i]
		V = []
		short_pulse = []
		spike = []
		syn = []
		noise1 = []
		noise2 = []

		while True: #split 1 picture to spyke,synapse,short_pulse and noise
			try:
				idx = np.where(t2 > -10)[0]
				if len(idx)==0:
					break
				idx=idx[0]+25000
				#V.append(t2[:idx])
				short_pulse.append(t2[69000:79000])
				syn.append(t2[130000:idx])
				spike.append(t2[110000:130000])
				#noise1.append(t2[:69000])
				noise2.append(t2[79000:110000])
				t2 = t2[idx:]

			except:
				break
		print('sp', np.shape(short_pulse),
			'syn', np.shape(syn),
			'spike', np.shape(spike),
			'noise2', np.shape(noise2),
			't2', np.shape(t2))

		syn = reshape_data(syn)
	#	V = reshape_data(V)
		with open(base+'syn/_'+str(index+i)+'.p', 'wb') as f:
			pickle.dump({'syn_'+str(index+i): np.array(syn)},f)
		with open(base+'short_pulse/_'+str(index+i)+'.p', 'wb') as f:
			pickle.dump({'short_pulse_'+str(index+i): np.array(short_pulse)}, f)
		with open(base + 'spike/_' + str(index + i) + '.p', 'wb') as f:
			pickle.dump({'spike_' + str(index + i): np.array(spike)}, f)
		with open(base + 'noise/_' + str(index + i) + '.p', 'wb') as f:
			pickle.dump({'noise2_' + str(index + i): np.array(noise2)}, f)
		del t2,t1,t,syn,spike,noise2,short_pulse,idx
		gc.collect()


		#	'noise1_'+str(index+i): np.array(noise1)),
		#}, open(base+'noise/1_'+str(index+i)+'.p', 'wb'))
		#pickle.dump({


	#return V,short_pulse,syn,spike,noise1,noise2
	#return short_pulse,syn,spike,noise2

