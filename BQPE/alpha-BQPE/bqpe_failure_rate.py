import random
import numpy as np
from numpy import pi, exp, cos
import os
import time as time
import sys
from datetime import datetime

from bqpe import *

M_range = np.linspace(1, 100, 100, dtype = int)
Attempts = 100

pres = 5*10**-3
MaxR = 1/pres**2

from progress.bar import ShadyBar
bar = ShadyBar('Generating:', max = Attempts * len(M_range), suffix = '%(percent).2f%%')
data = []

Alpha_List = -np.log(M_range)/np.log(pres)
Alpha_Values = [min(1, x) for x in Alpha_List]

for alpha in Alpha_Values:
    failures = 0
    idx = 0
    while idx < Attempts:
        r = random.uniform(-pi, pi)
        # start = time.time()
        flag, est, error, runs, sig = bqpe_analytical(threshold = pres, Phi = r, Alpha = alpha, sigma = pi / 4, Max_Runs = MaxR)
        # end = time.time()
        if flag == 1:
            failures+=1

        idx+=1
        bar.next()
        # else:
            # print('failed')
    data.append(failures/Attempts)
bar.finish()

# print(data)
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "helvetica"

plt.plot(M_range, data, linewidth = 2)
plt.grid(True)
plt.xlabel('Max M', fontsize = 15)
plt.ylabel('Rate of failure', fontsize = 15, labelpad = 5)
plt.xticks(fontsize = 12)
plt.yticks(fontsize = 12)
plt.legend(fontsize = 12)
plt.savefig('FailureRate_M_Uncapped.png', bbox_inches = 'tight')
plt.clf()