import random
import numpy as np
from numpy import pi, exp
import os
import time as time
import sys

def transpose(L):
    return(list(map(list,zip(*L))))

def aqpe(threshold = 5*(10**-3), Phi = 0, Alpha = 0, sigma = pi / 4, Sample_Size = 100, Max_Runs = 100000):
    if not 0<= Alpha<=1 or not -pi<=Phi<=pi:
        raise ValueError("Alpha: must be between 0 and 1 inclusive.")
    if sigma<0 or sigma > 1:
        raise ValueError("Sigma outside of acceptable ranges")

    mu = random.uniform(-pi,pi)
    run = 0
    flag = 0
    while sigma > threshold:
        Sampled = np.random.normal(mu, sigma, Sample_Size)
        
        M = max(1,int(round(1/(sigma**Alpha))))
        theta = mu - sigma
        p = 1/2 + np.cos(M*theta -  M*Phi)/2

        if random.uniform(0, 1) < p:
            outcome = 0
        else:
            outcome = 1
        
        accepted = []         
        for x_ in Sampled:
            P = 1/2 + (1-2*outcome)*np.cos(M*x_- M*theta)
            if P > random.uniform(0, 1):
                accepted.append(
                    x_
                )
        # print(len(accepted), P, p)
        if len(accepted) < 2:
            continue
        
        mu, sigma = np.mean(accepted), np.std(accepted)
        
        run += 1
        if run > Max_Runs:
            flag = 1
            break

    err = abs(
    abs(np.cos(Phi/2)) - abs(np.cos(mu/2))
    )
    return(flag,   mu, err, run, sigma)
    
r = random.uniform(-pi, pi)
print(
    aqpe(threshold = 5*(10**-3), Phi = r, Alpha = 0, sigma = pi / 4, Sample_Size = 1000, Max_Runs = 100000)
)