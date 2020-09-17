import random
import numpy as np
from numpy import pi, exp, cos, sin
import os
import time as time
import sys

def transpose(L):
    return(list(map(list,zip(*L))))

def Update_Prior(M, theta, outcome, mu, sigma):
    d = outcome
    
    Expectation = mu - (
        (1-2*d)*M*sigma**2 * sin(M*mu) * cos(M*theta)
        ) / (
            exp(M**2 * sigma**2 / 2) + (1-2*d)*cos(M*mu)*cos(M*theta)
            )

    VarNum = exp(M**2 * sigma**2) + (
        0.5*(2*d - 1)*cos(M*theta) * (
            2*exp(M**2 * sigma**2 / 2)*(M**2 * sigma**2 - 2)*cos(M*mu) + (2*d - 1)*cos(M*theta)*(
                1 - (2 * M**2 * sigma**2) + cos(2*M*mu)
            )
        )
        )

    VarDenom = (
        exp(M**2 * sigma**2 / 2) + (1 - 2*d)*cos(M*theta)*cos(M*mu)
        )**2

    Variance = sigma**2 * (VarNum / VarDenom)
    
    Std = np.sqrt(Variance)
    return(Expectation, Std)


def bqpe_analytical(threshold = 5*(10**-3), Phi = 0, MaxM = 1, sigma = pi / 4, Max_Runs = 10**5):
    if not -pi<=Phi<=pi:
        raise ValueError("Alpha: must be between 0 and 1 inclusive.")
    if sigma<0 or sigma > 1:
        raise ValueError("Sigma outside of acceptable ranges")

    mu = random.uniform(-pi, pi)
    run = 0
    flag = 0
    while sigma > threshold:
        
        M = max(
            1, min(
                np.floor(1/sigma + 1/2), MaxM
            )
        )
        p = 1/2 + cos(M*Phi)/2

        if random.uniform(0, 1) < p:
            outcome = 0
        else:
            outcome = 1
        
        mu, sigma = Update_Prior(M, 0, outcome, mu, sigma)
        
        run += 1
        if run > Max_Runs:
            flag = 1
            break

    err = abs(
    abs(cos(Phi/2)) - abs(cos(mu/2))
    )

    return(flag, float('%.5f'%(cos(mu/2))), err, run, sigma)



def bqpe_numerical(threshold = 5*(10**-3), Phi = 0, MaxM = 1, sigma = pi / 4, Sample_Size = 100, Max_Runs = 100000):
    if not -pi<=Phi<=pi:
        raise ValueError("Alpha: must be between 0 and 1 inclusive.")
    if sigma<0 or sigma > 1:
        raise ValueError("Sigma outside of acceptable ranges")

    mu = random.uniform(-pi, pi)
    sigma = sigma
    run = 0
    flag = 0
    while sigma > threshold:
        Sampled = np.random.normal(mu, sigma, Sample_Size)
        M = max(
            1, min(
                np.floor(1/sigma + 1/2), MaxM
            )
        )

        p = 1/2 + cos(M*Phi)/2

        if random.uniform(0, 1) < p:
            outcome = 0
        else:
            outcome = 1
        
        accepted = []         
        for varphi in Sampled:
            P = 1/2 + (1-2*outcome)*cos(M*varphi)/2
            if P > random.uniform(0, 1):
                accepted.append(
                    varphi
                )

        if len(accepted) < 2:
            sigma *= 1.2
            continue
        # mu, sigma = Update_Prior(M, theta, outcome, mu, sigma)
        mu, sigma = np.mean(accepted), np.std(accepted)
        
        run += 1
        if run > Max_Runs:
            flag = 1
            break

    err = abs(
    abs(cos(Phi/2)) - abs(cos(mu/2))
    )

    return(flag, float('%.5f'%(cos(mu/2))), err, run, sigma)

# for _ in range(50):

#     print(
#         bqpe_analytical(Phi = 0.324234234)
#     )