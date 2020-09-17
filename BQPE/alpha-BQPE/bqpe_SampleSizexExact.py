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
    raise ValueError('Needs input argument of 0, 1.')


L = 500
thres = 5*10**-3
# thres = 10**-3
MaxR = 10**5

current = datetime.now()


current = datetime.now()

analytical_write = 'data/'+prefix+'analytical_alpha_Tol%s.txt'%thres
numerical_write = 'data/'+prefix+'numerical_alpha_Tol%s.txt'%thres

write_here = prefix + 'analytical_vs_numerical_%s.txt'%current.strftime('%H%M-%d%m%y')
f = open(write_here,'w+')
f.write('Runs: %s\nTolerance: %s\n'%(L,thres))
f.close()

ana_time = []
ana_err = []
from progress.bar import ShadyBar
bar = ShadyBar('Generating analytical data... ', max = L, suffix = '%(percent).2f%%')
idx = 0
Failure_Rate = 0
while idx < L:
    r = random.uniform(-pi, pi)
    start = time.time()
    flag, est, error, runs, sig = bqpe_analytical(threshold = thres, Phi = r, Alpha = 0, sigma = pi / 4, Max_Runs = MaxR)
    end = time.time()
    if flag == 0:
        ana_time.append(end-start)
        ana_err.append(error)
        idx += 1
        bar.next()
    else:
        Failure_Rate += 1
        # print('Failed\n')
bar.finish()

f = open(write_here, "a+")
f.write('Analytical data:\n    Error = %.8f\n    Time = %.8fs\n'%(np.median(ana_err), np.median(ana_time)))
f.write('Failed runs: %s\n\n'%Failure_Rate)
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
bar = ShadyBar('Generating numerical data... ', max = L*len(sample_sizes), suffix = '%(percent).2f%%')
for S in sample_sizes:
    num_time = []
    num_err = []
    idx = 0
    Failure_Rate = 0
    while idx < L:
        r = random.uniform(-pi, pi)
        start = time.time()
        flag1, est1, error1, runs1, sig = bqpe_numerical(threshold = thres, Phi = r, Alpha = 0, sigma = pi / 4, Sample_Size = S, Max_Runs = MaxR)
        end = time.time()
        if flag1 == 0:
            num_time.append(end-start)
            num_err.append(error1)
            idx += 1
            bar.next()
        else:
            Failure_Rate += 1
    f = open(write_here, "a+")
    str1 = '  {:>11}'.format('%s'%S)
    str2 = '  {:>12}'.format('%.8f'%(np.median(num_err)))
    str3 = '{:>12}\n'.format('%.8f'%(np.median(num_time)))
    f.write(str1+'|'+str2+str3)
    
    f.close()
bar.finish()
f = open(write_here, "a+")
f.write('Failed runs: %s'%Failure_Rate)
f.close()