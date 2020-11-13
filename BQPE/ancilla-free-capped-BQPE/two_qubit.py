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

def circuit_no_ancilla(param = 0, sigma = 0, M = 1):  
    
    Engine = MainEngine()
    
    q1 = Engine.allocate_qureg(1)
    q2 = Engine.allocate_qureg(1)

    Rx(param) | q1
    Rx(param) | q2

    # [R * c-Pi * R^\dag * P * R * c-Pi * R^\dag * P]
    for _ in range(M):
        Z | q1
        Z | q2

        DaggeredGate(Rx(param)) | q1 
        DaggeredGate(Rx(param)) | q2
        
        Rzz(pi/2) | (q1, q2)

        Rz(pi/2) | q1
        Rz(pi/2) | q2

        Rx(param) | q1
        Rx(param) | q2

        Z | q1
        Z | q2

        DaggeredGate(Rx(param)) | q1 
        DaggeredGate(Rx(param)) | q2

        Rzz(pi/2) | (q1, q2)
        
        Rz(pi/2) | q1
        Rz(pi/2) | q2
   
        Rx(param) | q1
        Rx(param) | q2


    DaggeredGate(Rx(param)) | q1 
    DaggeredGate(Rx(param)) | q2

    Measure | q1
    Measure | q2

    Engine.flush()
    
    return(int(q1), int(q2))

m_range = np.linspace(1, 2, 2, dtype=int)
avg_over = 10000

Phi = pi/np.sqrt(2) #np.random.uniform(-pi, pi)
expectation = abs(cos(Phi / 2))
theta = np.arccos(np.sqrt(expectation))

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

# for m in m_range:
#     p_1 = 0
#     print('M = %s'%m)
#     ana_prob_1 = (1 - cos(2*m*Phi))/2
#     for _ in range(avg_over):
#         out1, out2 = circuit_no_ancilla(param = theta, M = m)
#         if out1 + out2 == 1:
#             p_1 += 1/avg_over
#     print("    Numerics: %.4f / %.4f"%(1-p_1, p_1) )
#     print("    Analytical: %.4f / %.4f"%(1-ana_prob_1, ana_prob_1))