import random
import numpy as np
from numpy import pi, exp, cos
import os
import time as time
import sys
from datetime import datetime
from bqpe import *



L = 10**3
thres = 5*10**-3
thres = 10**-3
max = 1/thres**2 + 0.1 * 1/thres

random_phases = []
for _ in range(L):
    random_phases.append(
        random.uniform(-pi, pi)
    )

current = datetime.now()



write_here = 'analytical_vs_numerical_%s.txt'%current.strftime('%H%M-%d%m%y')
f = open(write_here,'w+')
f.write('Runs: %s\nTolerance: %s\n'%(L,thres))
f.close()

ana_time = 0
ana_err = 0
for idx in range(L):
    r = random_phases[idx]
    start = time.time()
    flag, est, error, runs, sig = bqpe_analytical(threshold = thres, Phi = r, Alpha = 0, sigma = pi / 4, Max_Runs = max)
    end = time.time()
    
    ana_time += (end-start) / L
    ana_err += error / L

f = open(write_here, "a+")
f.write('Analytical data:\n    Error = %.8f\n    Time = %.8fs\n\n'%(ana_err, ana_time))
f.close()

f = open(write_here, "a+")
base_str = 'Numerical data:\n'
str_1 = '  {:<11}'.format('Sample Size')
str_2 = '  {:^12}'.format('Error')
str_3 = '{:^12}\n'.format('Time')
f.write(base_str+str_1+'|'+str_2+str_3)
f.write('{:_^40}\n'.format(''))
f.close()

sample_sizes = [10, 50, 100, 500, 1000]
for S in sample_sizes:
    num_time = 0
    num_err = 0
    for _ in range(L):
        start = time.time()
        flag1, est1, error1, runs1, sig = bqpe_numerical(threshold = thres, Phi = r, Alpha = 0, sigma = pi / 4, Sample_Size = S, Max_Runs = max)
        end = time.time()
        
        num_time += (end-start)/L
        num_err += error1/L

    f = open(write_here, "a+")
    str1 = '  {:>11}'.format('%s'%S)
    str2 = '  {:>12}'.format('%.8f'%num_err)
    str3 = '{:>12}\n'.format('%.8f'%num_time)
    f.write(str1+'|'+str2+str3)
    f.close()