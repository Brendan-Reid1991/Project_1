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
    phi_mu= mu+((1-2*d)*M*sigma**2*np.sin(M*mu - M*theta))/(np.exp(0.5*M**2*sigma**2)+(1-2*d)*np.cos(M*mu - M*theta))
    numa = (2*np.exp(M**2*sigma**2)+2*(2*d-1)*exp(M**2 * sigma**2 * 0.5)*(M**2*sigma**2-2)*np.cos(M*mu - M*theta)+(1-2*d)**2*(1-2*M**2*sigma**2+np.cos(2*M*mu - 2*M*theta)))
    demonetor=2*(np.exp(0.5*M**2*sigma**2)+(1-2*d)*np.cos(M*mu - M*theta))**2
    sigma_2=sigma**2*numa/demonetor
    phi_sigma=np.sqrt(sigma_2)
    
    return(phi_mu, phi_sigma)


def bqpe_analytical(threshold = 5*(10**-3), Phi = 0, Alpha = 0, sigma = pi / 4, Max_Runs = 10**5):
    if not 0<= Alpha<=1 or not -pi<=Phi<=pi:
        raise ValueError("Alpha: must be between 0 and 1 inclusive.")
    if sigma<0 or sigma > 1:
        raise ValueError("Sigma outside of acceptable ranges")

    mu = random.uniform(-pi, pi)
    run = 0
    flag = 0
    while sigma > threshold:
        theta = mu - sigma
        M = max(
            1, np.floor((1/sigma**Alpha) + 1/2)
        )

        p = 1/2 + cos(M*Phi - M*theta)/2

        if random.uniform(0, 1) < p:
            outcome = 0
        else:
            outcome = 1
        
        mu, sigma = Update_Prior(M, theta, outcome, mu, sigma)
        # mu, sigma = np.mean(accepted), np.std(accepted)
        
        run += 1
        if run > Max_Runs:
            flag = 1
            break

    err = abs(
    abs(cos(Phi/2)) - abs(cos(mu/2))
    )

    return(flag, float('%.5f'%(cos(mu/2))), err, run, sigma)



def bqpe_numerical(threshold = 5*(10**-3), Phi = 0, Alpha = 0, sigma = pi / 4, Sample_Size = 100, Max_Runs = 100000):
    if not 0<= Alpha<=1 or not -pi<=Phi<=pi:
        raise ValueError("Alpha: must be between 0 and 1 inclusive.")
    if sigma<0 or sigma > 1:
        raise ValueError("Sigma outside of acceptable ranges")

    mu = random.uniform(-pi, pi)
    sigma = sigma
    run = 0
    flag = 0
    while sigma > threshold:
        Sampled = np.random.normal(mu, sigma, Sample_Size)
        theta = mu - sigma
        M = max(
            1, np.floor((1/sigma**Alpha) + 1/2)
        )

        p = 1/2 + cos(M*Phi - M*theta)/2

        if random.uniform(0, 1) < p:
            outcome = 0
        else:
            outcome = 1
        
        accepted = []         
        for varphi in Sampled:
            P = 1/2 + (1-2*outcome)*cos(M*varphi-M*theta)/2
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