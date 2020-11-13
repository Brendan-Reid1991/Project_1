import random
import numpy as np
from numpy import pi, exp
import os
import time as time
import sys

def transpose(L):
    return(list(map(list,zip(*L))))





def Update_Prior(M,  measurements,mu,sigma):
    d= measurements
    phi_mu= mu-((1-2*d)*M*sigma**2*np.sin(M*mu))/ (np.exp(0.5*M**2*sigma**2)+(1-2*d)*np.cos(M*mu))
    numa = (np.exp(M**2*sigma**2))+0.5*(2*d-1)*(2*np.exp(0.5*M**2*sigma**2)*(M**2*sigma**2-2)*np.cos(M*mu)+(2*d-1)*(1-2*M**2*sigma**2+np.cos(2*M*mu)))
    demonetor=(np.exp(0.5*M**2*sigma**2)+(1-2*d)*np.cos(M*mu))**2
    sigma_2=sigma**2*numa/demonetor
    phi_sigma=np.sqrt(sigma_2)
    
    
    return(phi_mu, phi_sigma)



def aqpe(threshold = 5*(10**-3), D_max = 1,  mu = pi/2, sigma = pi / 4, Max_Runs = 10000, Phi = 0):

    run = 0
    while sigma > threshold:
        M = int(round(1/(2*sigma)))*2
        if M < 1:
            M = 1
        if M>D_max:
            M=D_max
        P = 1/2*(1+np.cos(M*Phi))
        x = random.uniform(0, 1)
        if x < P:
            measures = 0
        else: 
            measures = 1
        mu, sigma = Update_Prior(M, measures, mu, sigma)
        run += 1
        if run > Max_Runs:
            break
        #print(sigma, mu, M)

    err = abs(
    abs(np.cos(Phi/2)) - abs(np.cos(mu/2))
    )
    if run > Max_Runs:
        return(1, mu, err, run, M)
    else:
        return(0, mu, err, run, M)

from progress.bar import ShadyBar

Phi_Range = np.linspace(-pi, pi, 100)
average_over = 50
Results = []
bar = ShadyBar('Processing...', max = len(Phi_Range)*average_over, suffix = "%(percent).2f%%")

for phi in Phi_Range:
    ests = []
    for k in range(average_over):
        result=aqpe(threshold=5*10**-3, D_max=100,Phi=phi, sigma=pi/4, Max_Runs=10**5)
        ests.append(np.cos(result[1]/2))
        bar.next()
    Results.append(np.median(ests))
    # print("%.4f -- %.4f"%(phi, np.median(_errs)))
        # print(l,j, result[1], result[2],result[3], file=open("brendan_test.txt", "a"))
bar.finish()
import matplotlib.pyplot as plt
plt.plot(Phi_Range/pi, np.cos(np.asarray(Phi_Range)/2), linewidth = 1.5, color = "grey", label = "True")
plt.plot(Phi_Range/pi, Results, linewidth = 1.5, label = 'Concilla', linestyle = 'dotted')
plt.xlabel(r"$\phi/\pi$")
plt.ylabel(r"$\cos(\mu/2)$s")
plt.legend()
plt.grid(True)
plt.savefig('james_code.pdf')