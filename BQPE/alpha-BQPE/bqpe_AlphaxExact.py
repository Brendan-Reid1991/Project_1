import random
import numpy as np
from numpy import pi, exp, cos
import os
import time as time
import sys
from datetime import datetime

if sys.argv[1] == '0':
    from bqpe import *
    prefix = 'Superposition_'
elif sys.argv[1] == '1':
    from bqpe_collapsed import *
    prefix = 'Collapsed_'
else:
    raise ValueError('Needs integer input 0,1')


def transpose(L):
    return(list(map(list,zip(*L))))

L = 250
thres = 5*10**-3
# thres = 10**-3
MaxR = 10**5#2/thres**2 + 0.1 * 1/thres


current = datetime.now()

analytical_write = 'data/'+prefix+'analytical_alpha_Tol%s.txt'%thres
numerical_write = 'data/'+prefix+'numerical_alpha_Tol%s.txt'%thres

if os.path.exists(analytical_write):
    os.rename(analytical_write, 'data/analytical_alpha_Tol%s_Renamed%s.txt'%(thres, current.strftime('%d%m%y')))
if os.path.exists(numerical_write):
    os.rename(numerical_write, 'data/numerical_alpha_Tol%s_Renamed%s.txt'%(thres, current.strftime('%d%m%y')))

Alpha_Values = np.linspace(0, 1, 25)
from progress.bar import ShadyBar
bar = ShadyBar('Generating analytical data... ', max = L * len(Alpha_Values), suffix = '%(percent).2f%%')
analytical_data = []
for alpha in Alpha_Values:
    ana_time = []
    ana_err = []
    idx = 0
    while idx < L:
        r = random.uniform(-pi, pi)
        start = time.time()
        flag, est, error, runs, sig = bqpe_analytical(threshold = thres, Phi = r, Alpha = alpha, sigma = pi / 4, Max_Runs = MaxR)
        end = time.time()
        if flag == 0:
            ana_time.append(end-start)
            ana_err.append(error)
            idx += 1
            bar.next()
        # else:
            # print('failed')
    analytical_data.append([alpha, np.median(ana_err), np.median(ana_time)])
bar.finish()
np.savetxt(analytical_write, analytical_data)

bar = ShadyBar('Generating numerical data... ', max = L * len(Alpha_Values), suffix = '%(percent).2f%%')
numerical_data = []
for alpha in Alpha_Values:
    num_time = []
    num_err = []
    idx = 0
    while idx < L:
        r = random.uniform(-pi, pi)
        start = time.time()
        flag, est, error, runs, sig = bqpe_numerical(threshold = thres, Phi = r, Alpha = alpha, sigma = pi / 4, Max_Runs = MaxR, Sample_Size = 500)
        end = time.time()
        if flag == 0:
            num_time.append(end-start)
            num_err.append(error)
            idx += 1
            bar.next()
    numerical_data.append([alpha, np.median(num_err), np.median(num_time)])
bar.finish()
np.savetxt(numerical_write, numerical_data)


import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "helvetica"

plt.plot(Alpha_Values, transpose(analytical_data)[1], linewidth = 2, label = 'Analytical')
plt.plot(Alpha_Values, transpose(numerical_data)[1], linewidth = 2, label = 'Numerical')
plt.grid(True)
plt.xlabel(r'$\alpha$', fontsize = 15)
plt.ylabel('Median Error', fontsize = 15, labelpad = 5)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.legend(fontsize = 12)
plt.savefig(prefix+'Analytical_Numerical_ErrorAgainstAlpha.png', bbox_inches = 'tight')
plt.clf()

plt.plot(Alpha_Values, transpose(analytical_data)[2], linewidth = 2, label = 'Analytical')
plt.plot(Alpha_Values, transpose(numerical_data)[2], linewidth = 2, label = 'Numerical')
plt.grid(True)
plt.xlabel(r'$\alpha$', fontsize = 15)
plt.ylabel('Median Timing', fontsize = 15, labelpad = 5)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.legend(fontsize = 12)
plt.savefig(prefix + 'Analytical_Numerical_TimeAgainstAlpha.png', bbox_inches = 'tight')
plt.clf()