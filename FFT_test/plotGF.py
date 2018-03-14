import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import re
from scipy.stats import norm
import os
import glob
import re
import linecache

mpl.style.use('seaborn')
sns.set_style("whitegrid")

regD = r"([0-9]+)\_"
regN = r"\_(.+)\_"

data_tails = []
data_mf = []
data_mf_names = []
data_it = []
data_it_names = []

#MF: 
#	 mFreq 	Spin UP Re 	 Spin UP Im 	 Spin DOWN Re 	 Spin DOWN Im 	Spin UP Re Err 	 Spin UP Im Err 	 Spin DOWN Re Err 	 Spin DOWN Im Err 
#it
#iTime 	Spin up 	Spin down 	Spin up Err	Spin down Err
def argsortll(seq, index):
    # http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
    return sorted(range(len(seq)), key=seq.__getitem__(index))


for filename in glob.glob(r'*.out'):
    ii = 0
    if filename[-6:-4] == 'MF':
        ii = 0
    elif filename[-6:-4] == 'IT':
        ii = 1
    else:
        continue;
    with open(filename, 'r') as infile:
        it = int(re.search(regD, filename).group(0)[:-1])
        t = re.search(regN, filename).group(0)[1:-1]
        if ii == 0:
            dat_tail = []
            line = infile.readline()
            while line.startswith("#") or len(dat_tail) < 2:
                if not line.startswith("#"):
                    dat_tail.append(np.fromstring(line, sep=' '))
                line = infile.readline()
            dat = np.genfromtxt(infile, skip_header = 7)
            data_mf_names.append([t,it])
            data_mf.append(dat)
            data_tails.append(dat_tail)
        if ii == 1:
            dat = np.genfromtxt(infile, skip_header = 1)
            data_it_names.append([t,it])
            data_it.append(dat)

data_tails = np.array(data_tails)
data_mf = np.array(data_mf)
data_it = np.array(data_it)
data_mf_names = np.array(data_mf_names)
data_it_names = np.array(data_it_names)
#i_mf = data_mf_names[:,1].argsort()
#i_it = data_it_names[:,1].argsort()
#data_tails = data_tails[i_mf]
#data_mf = data_mf[i_mf]
#data_it = data_it[i_it]
#data_mf_names = data_mf_names[i_mf]
#data_it_names = data_it_names[i_it]


err_mf = np.zeros((data_mf_names.shape[0], data_mf.shape[1]))
err_it = np.zeros((data_mf_names.shape[0], data_mf.shape[1]))
for i, it in enumerate(data_mf_names[:,1]):

    err_mf[it] = 0
    err_it[it] = 0
