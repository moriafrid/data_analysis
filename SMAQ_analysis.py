from open_pickle import read_from_pickle
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

EPSP1=read_from_pickle('/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/data/syn/syn.p')

EPSP=read_from_pickle('/ems/elsc-labs/segev-i/moria.fridman/project/data_analysis_git/data_analysis/data/correct_syn/clear_syn.p')
higth=np.max(np.array(EPSP[0]),axis=1)
n, bins_hist, patches = plt.hist(higth, bins=12, color='#0504aa',
                            alpha=0.7, rwidth=0.85)
plt.show()

mean = (higth+77.3).mean()
var = (higth+77.3).var()
skew = stats.skew(higth+77.3)
std = (higth+77.3).std()
N = -mean**2/((skew*mean-std)*std)
counts, bins = np.histogram(higth)

mids = 0.5*(bins[1:] + bins[:-1])
probs = counts / np.sum(counts)

mean = np.sum(probs * mids)
sd = np.sqrt(np.sum(probs * (mids - mean)**2))
