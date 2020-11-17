import random
import numpy as np
from numpy import pi, exp, sin, cos
from projectq import MainEngine
import os
import time as time
import sys
from projectq.ops import (H, R, Rz, Rx, Ry, C, Measure, DaggeredGate, Z, X, Y, Rzz, MatrixGate, All)

def circuit_no_ancilla(param = 0, M = 1):  
    
    Engine = MainEngine()
    
    q1 = Engine.allocate_qureg(1)

    Rx(param) | q1
    # M = max(1, int(M/2))
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

def p0(phi = 0, M = 1):
    return(
        1/2 + cos(2*M*phi)/2
    )

Exp = 0.3453454
theta = np.arccos(Exp)
Phi = 2*theta
avg = 500

for m in range(1, 11):
    numerical_p1 = 0
    for _ in range(5000):
        out = circuit_no_ancilla(param = theta, M = m)
        numerical_p1 += out/5000
    print(m, "%.4f || %.4f"%(1 - numerical_p1, p0(phi = Phi, M = m)))