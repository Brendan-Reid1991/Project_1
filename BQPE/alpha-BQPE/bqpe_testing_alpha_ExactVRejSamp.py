import random
import numpy as np
from numpy import pi, exp, cos
import os
import time as time
import sys
from datetime import datetime
from bqpe import *

def transpose(L):
    return(list(map(list,zip(*L))))

L = 500
thres = 5*10**-3
# thres = 10**-3
MaxR = 10**5#2/thres**2 + 0.1 * 1/thres


random_phases = []
for _ in range(L):
    random_phases.append(
        random.uniform(-pi, pi)
    )

current = datetime.now()

analytical_write = 'data/analytical_alpha_Tol%s.txt'%thres
numerical_write = 'data/numerical_alpha_Tol%s.txt'%thres

if os.path.exists(analytical_write):
    os.rename(analytical_write, 'data/analytical_alpha_Tol%s_Renamed%s.txt'%(thres, current.strftime('%d%m%y')))
if os.path.exists(numerical_write):
    os.rename(numerical_write, 'data/numerical_alpha_Tol%s_Renamed%s.txt'%(thres, current.strftime('%d%m%y')))

Alpha_Values = np.linspace(0, 1, 25)
from progress.bar import ShadyBar
bar = ShadyBar('Generating analytical data... ', max = L * len(Alpha_Values), suffix = '%(percent).2f%%')
analytical_data = []
for alpha in Alpha_Values:
    ana_time = 0
    ana_err = 0
    idx = 0
    while idx < L:
        r = random_phases[idx]
        start = time.time()
        flag, est, error, runs, sig = bqpe_analytical(threshold = thres, Phi = r, Alpha = alpha, sigma = pi / 4, Max_Runs = MaxR)
        end = time.time()
        if flag == 0:
            ana_time += (end-start) / L
            ana_err += error / L
            idx += 1
            bar.next()
        # else:
            # print('failed')
    analytical_data.append([alpha, ana_err, ana_time])
bar.finish()
np.savetxt(analytical_write, analytical_data)

bar = ShadyBar('Generating numerical data... ', max = L * len(Alpha_Values), suffix = '%(percent).2f%%')
numerical_data = []
for alpha in Alpha_Values:
    num_time = 0
    num_err = 0
    idx = 0
    while idx < L:
        r = random_phases[idx]
        start = time.time()
        flag, est, error, runs, sig = bqpe_numerical(threshold = thres, Phi = r, Alpha = alpha, sigma = pi / 4, Max_Runs = MaxR, Sample_Size = 500)
        end = time.time()
        if flag == 0:
            num_time += (end-start) / L
            num_err += error / L
            idx += 1
            bar.next()
    numerical_data.append([alpha, num_err, num_time])
bar.finish()
np.savetxt(numerical_write, numerical_data)


import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "helvetica"

plt.plot(Alpha_Values, transpose(analytical_data)[1], linewidth = 2, label = 'Analytical')
plt.plot(Alpha_Values, transpose(numerical_data)[1], linewidth = 2, label = 'Numerical')
plt.grid(True)
plt.xlabel(r'$\alpha$', fontsize = 15)
plt.ylabel('Average Error', fontsize = 15, labelpad = 5)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.legend(fontsize = 12)
plt.savefig('Analytical_Numerical_ErrorAgainstAlpha.png', bbox_inches = 'tight')
plt.clf()

plt.plot(Alpha_Values, transpose(analytical_data)[2], linewidth = 2, label = 'Analytical')
plt.plot(Alpha_Values, transpose(numerical_data)[2], linewidth = 2, label = 'Numerical')
plt.grid(True)
plt.xlabel(r'$\alpha$', fontsize = 15)
plt.ylabel('Average Timing', fontsize = 15, labelpad = 5)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.legend(fontsize = 12)
plt.savefig('Analytical_Numerical_TimeAgainstAlpha.png', bbox_inches = 'tight')
plt.clf()