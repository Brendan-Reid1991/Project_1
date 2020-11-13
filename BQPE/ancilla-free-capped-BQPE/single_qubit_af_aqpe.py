import random
import numpy as np
from numpy import pi, exp, sin, cos
from projectq import MainEngine
import os
import time as time
import sys
from projectq.ops import (H, R, Rz, Rx, Ry, C, Measure, DaggeredGate, Z, X, Y, Rzz, MatrixGate, All)

def transpose(L):
    return(list(map(list,zip(*L))))

def ExactUpdate(outcome = 0, sigma = 0, mu = 0, M = 0):
    d = int(outcome)

    Expectation = mu - (
        (1-2*d)*M*sigma**2 * sin(M*mu)
        ) / (
            exp(M**2 * sigma**2 / 2) + (1-2*d)*cos(M*mu)
            )

    VarNum = exp(M**2 * sigma**2) + (
        0.5*(2*d - 1) * (
            2*exp(M**2 * sigma**2 / 2)*(M**2 * sigma**2 - 2)*cos(M*mu) + (2*d - 1)*(
                1 - (2 * M**2 * sigma**2) + cos(2*M*mu)
            )
        )
        )

    VarDenom = (
        exp(M**2 * sigma**2 / 2) + (1 - 2*d)*cos(M*mu)
        )**2

    Variance = sigma**2 * (VarNum / VarDenom)
    Std = np.sqrt(Variance)
    return(Expectation, Std)

Pauli = [Z]

def circuit_no_ancilla(param = 0, M = 1):  
    
    Engine = MainEngine()
    
    q1 = Engine.allocate_qureg(1)

    Rx(param) | q1

    # [R * c-Pi * R^\dag * P * R * c-Pi * R^\dag * P]
    for _ in range(M):
        Z | q1

        DaggeredGate(Rx(param)) | q1 

        Z | q1

        Rx(param) | q1


        Z | q1

        DaggeredGate(Rx(param)) | q1 

        Z | q1
   
        Rx(param) | q1


    DaggeredGate(Rx(param)) | q1 


    Measure | q1


    Engine.flush()
    
    return(int(q1))

m_range = np.linspace(1, 2, 2, dtype=int)
avg_over = 10000



Exp = 0.3453454
theta = np.arccos(Exp)
Phi = 2*theta
avg = 500
M = 10
for mm in range(1, 10):
    print("M = ", mm)
    analytical_p0 = 1/2*(1+cos(2*mm*Phi))
    num_p1 = 0
    for _ in range(avg):
        out = circuit_no_ancilla(param = theta, M = mm)
        num_p1 += out/avg
    print("    ",1-num_p1)
    print("    ",analytical_p0)
exit()


Expectations = []
from progress.bar import ShadyBar
bar = ShadyBar("Processing...", max = 100, suffix = "%(percent)d%%")
i = 0
while i < 100:
    mu = random.uniform(-pi, pi)
    sigma = pi/4
    run = 0
    while sigma > 0.005:
        M = int(np.round(1/sigma))
        out = circuit_no_ancilla(param = theta, M = M)
        mu, sig = ExactUpdate(outcome = out, sigma = sigma, mu = mu, M = 2*M)
        sigma = sig
        run += 1
        if run > 10**4:
            break
    if run > 10**4:
        pass
    else:
        Expectations.append(
            cos(mu/2)
        )
        i+=1
        bar.next()
print(np.median(Expectations))
print(np.mean(Expectations))
print(Expectations)



# for m in m_range:
#     counts = [0,0,0,0]
#     print('M = %s'%m)
#     for _ in range(avg_over):
#         out1, out2 = circuit_no_ancilla(param = theta, M = m)
#         p_1 = (1 - cos(2*m*Phi))/2
#         if out1 == 1:
#             if out2 == 1:
#                 counts[3] += 1/avg_over
#             else:
#                 counts[2] += 1/avg_over
#         else:
#             if out2 == 0:
#                 counts[0] += 1/avg_over
#             else:
#                 counts[1] += 1/avg_over
#     print("    00: %.4f | else: %.4f (01: %.4f, 10: %.4f, 11: %.4f)"%(counts[0], sum(counts[1:]),counts[1],counts[2],counts[3]))
#     print("    Analytical p0 = %.4f | p1 = %.4f"%(1-p_1, p_1))
