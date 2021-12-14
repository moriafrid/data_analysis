import matplotlib.pyplot as plt
import numpy as np
from open_pickle import read_from_pickle
import os
from add_figure import add_figure
import pickle

path='data/correct_syn'
try:os.mkdir(path)
except FileExistsError:pass
path1=path+'/syn2clear'
try:os.mkdir(path1)
except FileExistsError:pass

V_units,t_units=read_from_pickle('/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/data/syn/syn.p')
V=np.array(V_units)
temp_syn=np.mean(V,axis=0)
syn_time2clear1=np.argmax(temp_syn)-100
syn_time2clear2=np.argmax(temp_syn)+200
rest,new_syn=[],[]
for v in V:
    rest.append(np.mean(v[syn_time2clear1 - 100:syn_time2clear1]))
    v = v - np.mean(v[syn_time2clear1 - 100:syn_time2clear1])
    new_syn.append(v[syn_time2clear1 - 500:syn_time2clear2 + 1000])
REST=np.mean(rest)
t=t_units[syn_time2clear1 - 500:syn_time2clear2 + 1000]
for i,bolt_trace in enumerate(new_syn):
    add_figure('trace num '+str(i)+'\nmean on 100 points',str(syn_time2clear1 - 500)+':'+str(syn_time2clear2 + 1000),'mv')
    for v in new_syn:
        plt.plot(v,'grey', alpha=0.05,lw=0.5)
    mean_syn=np.mean(new_syn,axis=0)
    plt.plot(mean_syn,'black',lw=2)
    plt.plot(bolt_trace,'green',alpha=0.5,lw=1)
    print(i)
    plt.savefig(path1+'/trace_num'+str(i))
    plt.close()
    a=1


not_sure=[12,32,34,37]
# rigth=[0,2,3,4,5,6,7,8,9,10,14,15,17,18,20,21,22,23,24,25,27,28,33,35,36,40,41,45,47,48,49,50,51,53,54,55,56,57,58,59,62,65,66,68,70,71,72,73,74,76,78,79]#for time2syn+1000
# rigth=[0,2,3,4,5,6,7,8,9,10,12,14,15,17,18,20,21,22,23,24,25,27,28,32,33,35,36,40,41,45,47,48,49,50,51,53,54,55,56,57,58,59,62,65,66,68,70,71,72,73,74,76,78,79] #for time2syn+1750:
# false=[1,11,12,13,16,19,26,29,30,31,32,34,37,38,39,42,43,44,46,52,60,61,63,64,67,68,69,75,77,80,81,82,83,84,85,86,87,88,89,90,91]#for time2syn+1000

#for time2syn+1000 and stable antil point 1200
not_sure=[12,20,26,28,47,48,49,70,71,74,78,81,85,88]
cut_on_1000=[8,19,21,34,35,38,40,47,48,77,81,90]
rigth=[0,1,2,3,4,5,6,9,10,14,18,21,22,23,24,25,26,27,33,35,36,37,41,43,44,45,47,48,50,51,53,55,56,58,59,62,65,66,67,72,73,74,75,76,78,79,82,86] #for time2syn+1000 and stable antil point 1200
false=[7,11,13,15,16,17,29,39,42,46,49,52,54,57,60,61,63,64,68,69,71,80,83,84,85,87,88,89,91]#for time2syn+1000 and stable antil point 1200
path2=path+'/syn2clear_again'
try:os.mkdir(path2)
except FileExistsError:pass
new_syn1=[]
for num in rigth:
    new_syn1.append(new_syn[num])

for i,bolt_trace in enumerate(new_syn1):
    add_figure('trace num '+str(i)+'\nmean on 100 points',str(syn_time2clear1 - 500)+':'+str(syn_time2clear2 + 1000),'mv')
    for v in new_syn1:
        plt.plot(v,alpha=0.5)
    mean_syn=np.mean(new_syn1,axis=0)
    plt.plot(mean_syn,'black',linewidth=3)
    plt.plot(bolt_trace,'b',linewidth=2)
    print(i)
    plt.savefig(path2+'/trace_num'+str(i))
    plt.close()
    a=1
# remove=[18,23,27,34,43,47] #6? 14?

#for time2syn+1000 and stable antil point 1200
remove=[11,19,24,]

new_syn2=np.delete(new_syn1,remove,axis=0)+REST
add_figure('correct synapse',V_units.units,t_units.units)
for v in new_syn2:
    plt.plot(t,v, alpha=0.1,lw=0.5)
plt.plot(t,np.mean(new_syn2,axis=0),'black',lw=2)
plt.savefig(path+'/clear_syn')
with open(path + '/clear_syn.p', 'wb') as f:
    pickle.dump([new_syn2*V_units.units,t] , f)
add_figure('mean correct synepses',V_units.units,t_units.units)
plt.plot(t,np.mean(new_syn2,axis=0),'black',linewidth=2)
plt.savefig(path+'/mean_clear_syn')
with open(path + '/mean_syn.p', 'wb') as f:
    pickle.dump({'mean':[(np.mean(new_syn2, axis=0)) * V_units.units,t],'E_pas':REST}, f)

a=0
