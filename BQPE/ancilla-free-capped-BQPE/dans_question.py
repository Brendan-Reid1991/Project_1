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


