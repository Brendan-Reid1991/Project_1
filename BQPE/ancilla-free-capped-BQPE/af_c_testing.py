import random
import numpy as np
from numpy import pi, exp, sin, cos
from projectq import MainEngine
import os
import time as time

from af_capped_bqpe import *

M_max = 10**3
avg = 50
Estimate_Results = []
Error_Results = []
Run_Results = []
Phi_range = np.linspace(-pi, pi, 100)
for phi_val in Phi_range:
    i = 0
    median_est = []
    median_err = []
    median_run = []
    while i < avg:
        flag, exp, err, run, sigma = bqpe_analytical(threshold = 5*(10**-3), Phi = phi_val, MaxM = M_max, sigma = pi / 4, Max_Runs = 10**4)
        if flag == 1:
            pass
        else:
            median_est.append(exp)
            median_err.append(err)
            median_run.append(run)
            i+=1
    Estimate_Results.append(
        np.median(median_est)
    )
    Error_Results.append(
        np.median(median_err)
    )
    Run_Results.append(
        np.median(median_run)
    )

np.savetxt(np.asarray(Estimate_Results), "Estimate_AgainstPhi.txt")
np.savetxt(np.asarray(Error_Results), "Error_AgainstPhi.txt")
np.savetxt(np.asarray(Run_Results), "Run_AgainstPhi.txt")

import matplotlib as pyplot

plt.plot