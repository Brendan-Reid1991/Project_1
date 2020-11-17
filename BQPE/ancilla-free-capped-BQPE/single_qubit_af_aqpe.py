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
    for _ in range(int(M/2)):
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




Expectation_Value = random.uniform(-1, 1)  # Randomly chosen!
theta = np.arccos(Expectation_Value)  # Ansatz Parameter for this case
Phi = 2*theta # Corresponding Phase Value
avg = 500




print("\nCorrect statistics are replicated by the circuit!")
for M in range(1,11):
    # m=M     
    # if M%2 == 0 and M > 1:
    #     pass
    # else:
    #     M = max(1, M - 1)
    
    nump1 = 0
    for _ in range(10**3):
        out = circuit_no_ancilla(param = theta, M = 2*M)
        nump1 += out/10**3
    nump0 = 1-nump1
    P0 = 1/2 + cos(2*M*Phi)/2
    print(M,"||", P0,"|", nump0)
exit()




def circ_vqe(param = 0):
    Engine = MainEngine()
    
    q1 = Engine.allocate_qureg(1)

    Rx(param) | q1

    Measure | q1

    return(int(q1))


# p1 = 0
# for _ in range(avg):
#     p1 += circ_vqe(param = theta)
# p0 = 1-p1/avg

# print("Vqe",2*p0-1)

# print(p0)
# print(1/2*(1+cos(p0)))
# exit()

outcomes = []
i = 0
p1 = 0
Initial_Exp = 500
while i < Initial_Exp:
    out = circ_vqe(param = theta)
    p1 += out
    mu = 2*np.arccos(1 - 2*p1/(i + 1))
    outcomes.append(mu) 
    i+=1

mu = outcomes[-1]
sigma = np.std(outcomes)

print("After initial 500 shots of VQE, mu = %.4f and sigma = %.4f\nTrue value of mu should be %.4f"%(mu, sigma, Phi))
run = 0
while sigma > 0.001:
    M = max(1, round(0.5/sigma))
    if M%2 == 0 and M > 1:
        pass
    else:
        M = max(1, M - 1)

    out = circuit_no_ancilla(param = theta, M = M)

    mu, sig = ExactUpdate(outcome = out, sigma = sigma, mu = mu, M = M)
    sigma = sig
    run +=1
    if run>10**4:
        break

print("                   True  /  Estimated\nPhase phi::      %.5f / %.5f"%(Phi, mu))
print("Expectation::    %.5f / %.5f"%(cos(Phi / 2), cos(mu/2)))    
exit()

