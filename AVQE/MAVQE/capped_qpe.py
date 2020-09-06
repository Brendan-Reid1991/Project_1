import random
import numpy as np
from numpy import pi, exp
import os
import time as time
import sys



def Generate_Shots(Phi, M, theta, shot_number, Collapsed = 0):

    if Collapsed == 0:
        p = 1/2 + np.cos(M*theta)*np.cos(M*Phi)/2
    else:
        p = 1/2 + np.cos(M*theta-M*Phi)/2
    
    measurements = []
    
    for j in range(shot_number):
        if random.uniform(0, 1) < p:
            measurements.append(0)
        else:
            measurements.append(1)
    
    return(measurements)



def Probability(E, M, Theta, Phi, Collapsed = 0):
    if Collapsed == 0:
        return(
        1/2 + ((-1)**E * np.cos(M*Theta)*np.cos(M*Phi))/2
        )
    else:
        return(
        1/2 + ((-1)**E *np.cos(M*Theta-M*Phi))/2
        )


def Update_Prior(M, theta, measurement_outcomes, sampled, C = 0):

    N = int(len(measurement_outcomes))
    N_0 = int(measurement_outcomes.count(0))
    N_1 = N - N_0

    accepted = []
    recov = 1
    if N_0 == 0 or N_1 == 0:
        f0 = 1
        f1 = 1
    else:
        f0 = N/N_0
        f1 = N/N_1
    for phi in sampled:
        p0 = Probability(0, M, theta, phi, Collapsed = C)*f0
        p1 = Probability(1, M, theta, phi, Collapsed = C)*f1

        p = (p0**N_0)*(p1**N_1)

        if p >= random.uniform(0,1):
            accepted.append(phi)

    
    if len(accepted) < 2:
        phi_mu, phi_sigma = np.mean(sampled), (1+recov)*np.std(sampled)
    else:
        phi_mu, phi_sigma = np.mean(accepted), np.std(accepted)

    return(phi_mu, phi_sigma)



def aqpe_dynamic(Phi = 0, threshold = 5*(10**-3), Alpha = 0, nSamples = 100, nShots = 10, sigma = 3*pi / 4, Max_Runs = 1000, Col = 0):

    if Col not in [0, 1]:
        raise ValueError("Col: must be either 0 (superposition) or 1 (collapsed).")
    if not 0<= Alpha<=1:
        raise ValueError("Alpha: must be between 0 and 1 inclusive.")
    mu = np.random.uniform(-pi, pi)
    run = 0
    while sigma > threshold:
        M = int(round(1/(sigma**Alpha)))
        if M < 1:
            M = 1
        theta = mu - sigma
        samples = np.random.normal(mu, sigma, nSamples)
        measures = Generate_Shots(Phi, M, theta, nShots, Collapsed = Col)
        mu, sigma = Update_Prior(M, theta, measures, samples, C = Col)
        
        run += 1

        err = abs(
        abs(np.cos(Phi/2)) - abs(np.cos(mu/2))
        )

        if run > Max_Runs:
            break
        sys.stdout.write('\r %s -- %.4f  |  %.4f'%(nShots*run, sigma, err))
        sys.stdout.flush()
        # time.sleep(0.05)

    err = abs(
    abs(np.cos(Phi/2)) - abs(np.cos(mu/2))
    )
    if run > Max_Runs:
        return(1, err, nShots*run, M)
    else:
        return(0, err, nShots*run, M)



def aqpe(Phi = 0, threshold = 5*(10**-3), Alpha = 0, nSamples = 100, nShots = 10, sigma = 3*pi / 4, Max_Runs = 1000, Col = 0):

    if Col not in [0, 1]:
        raise ValueError("Col: must be either 0 (superposition) or 1 (collapsed).")
    if not 0<= Alpha<=1:
        raise ValueError("Alpha: must be between 0 and 1 inclusive.")
    mu = np.random.uniform(-pi, pi)
    run = 0
    while sigma > threshold:
        M = int(round(1/(sigma**Alpha)))
        if M < 1:
            M = 1
        theta = mu - sigma
        samples = np.random.normal(mu, sigma, nSamples)
        measures = Generate_Shots(Phi, M, theta, nShots, Collapsed = Col)
        mu, sigma = Update_Prior(M, theta, measures, samples, C = Col)
        run += 1
        if run > Max_Runs:
            break

    err = abs(
    abs(np.cos(Phi/2)) - abs(np.cos(mu/2))
    )
    if run > Max_Runs:
        return(1, err, nShots*run, M)
    else:
        return(0, err, nShots*run, M)

