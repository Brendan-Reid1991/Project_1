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


def Ansatz(param):
    return([
            [
                [Rx(param)]
            ]
            # ,[
            #     [Rx(param)]
            # ]
        ])


def AnsatzConj(param):
    Prep = Ansatz(param)
    PrepDag = []
    for sublist in Prep:
        new_sublist = []
        for subsublist in sublist:
            if len(subsublist) == 1:
                new_sublist.append([
                    DaggeredGate(subsublist[0])
                ])
            elif len(subsublist) == 2:
                new_sublist.append([
                    DaggeredGate(subsublist[0]), subsublist[1]
                ])

        PrepDag.append(
            new_sublist[::-1]
        )
    return(
        PrepDag
    )

Pauli = [Z]

def circuit(param = 0, sigma = 0, maxM = 10):  

    M = min(
        max(1, round(1/sigma)), maxM
    )

    Instruction_Set = Ansatz(param)
    Instruction_Set_Conjugated = AnsatzConj(param)


    n_qubits = len(Instruction_Set)
    
    Engine = MainEngine()

    mat = np.diag([1]*2*n_qubits)
    mat[0, 0] = mat[0, 0] - 2 
    Grover_Reflection = MatrixGate(mat)
    
    control = Engine.allocate_qubit()
    qubits = Engine.allocate_qureg(n_qubits)

    # Out-of-phase |+> state in control qubit
    H | control
    # R(- M * theta ) | control

    # State Preparation
    for idx,individual_instructions in enumerate(Instruction_Set):
        for gate_command in individual_instructions:
            L = len(gate_command)
            if L == 1:
                gate_command[0] | qubits[idx]
            elif L == 2:
                C(gate_command[0]) | (qubits[gate_command[1]] , qubits[idx])
            else:
                raise Exception('Instruction set on qubit %s has incorrect list length.'%idx)

    # Controlled - U implementation
    # [R * c-Pi * R^\dag * P * R * c-Pi * R^\dag * P]
    for _ in range(M):
        for idx, P in enumerate(Pauli):
            P | qubits[idx]

        # R^\dag
        for idx,individual_instructions in enumerate(Instruction_Set_Conjugated):
            for gate_command in individual_instructions:
                L = len(gate_command)
                if L == 1:
                    gate_command[0] | qubits[idx]
                elif L == 2:
                    C(gate_command[0]) | (qubits[gate_command[1]] , qubits[idx])
        # c-Pi
        C(Grover_Reflection) | (control, qubits)

        # R
        for idx,individual_instructions in enumerate(Instruction_Set):
            for gate_command in individual_instructions:
                L = len(gate_command)
                if L == 1:
                    gate_command[0] | qubits[idx]
                elif L == 2:
                    C(gate_command[0]) | (qubits[gate_command[1]] , qubits[idx])
        
        for idx, P in enumerate(Pauli):
            P | qubits[idx]

    # R^\dag
        for idx,individual_instructions in enumerate(Instruction_Set_Conjugated):
            for gate_command in individual_instructions:
                L = len(gate_command)
                if L == 1:
                    gate_command[0] | qubits[idx]
                elif L == 2:
                    C(gate_command[0]) | (qubits[gate_command[1]] , qubits[idx])
        # c-Pi
        C(Grover_Reflection) | (control, qubits)

    # R
        
        for idx,individual_instructions in enumerate(Instruction_Set):
            for gate_command in individual_instructions:
                L = len(gate_command)
                if L == 1:
                    gate_command[0] | qubits[idx]
                elif L == 2:
                    C(gate_command[0]) | (qubits[gate_command[1]] , qubits[idx])

    H | control

    Measure | control
    All(Measure) | qubits

    Engine.flush()

    return(
        int(control)
    )


def circuit_no_ancilla(param = 0, sigma = 0, maxM = 10):  

    M = min(
        max(1, round(1/sigma)), maxM
    )

    Instruction_Set = Ansatz(param)
    Instruction_Set_Conjugated = AnsatzConj(param)


    n_qubits = len(Instruction_Set)
    
    Engine = MainEngine()

    mat = np.diag([1]*2*n_qubits)
    mat[0, 0] = mat[0, 0] - 2 
    Grover_Reflection = MatrixGate(mat)
    
    qubits = Engine.allocate_qureg(n_qubits)


    # State Preparation
    for idx,individual_instructions in enumerate(Instruction_Set):
        for gate_command in individual_instructions:
            L = len(gate_command)
            if L == 1:
                gate_command[0] | qubits[idx]
            elif L == 2:
                C(gate_command[0]) | (qubits[gate_command[1]] , qubits[idx])
            else:
                raise Exception('Instruction set on qubit %s has incorrect list length.'%idx)

    # Controlled - U implementation
    # [R * c-Pi * R^\dag * P * R * c-Pi * R^\dag * P]
    for _ in range(M):
        for idx, P in enumerate(Pauli):
            P | qubits[idx]

        # R^\dag
        for idx,individual_instructions in enumerate(Instruction_Set_Conjugated):
            for gate_command in individual_instructions:
                L = len(gate_command)
                if L == 1:
                    gate_command[0] | qubits[idx]
                elif L == 2:
                    C(gate_command[0]) | (qubits[gate_command[1]] , qubits[idx])
        # c-Pi
        Grover_Reflection | qubits

        # R
        for idx,individual_instructions in enumerate(Instruction_Set):
            for gate_command in individual_instructions:
                L = len(gate_command)
                if L == 1:
                    gate_command[0] | qubits[idx]
                elif L == 2:
                    C(gate_command[0]) | (qubits[gate_command[1]] , qubits[idx])
        
        for idx, P in enumerate(Pauli):
            P | qubits[idx]

    # R^\dag
        for idx,individual_instructions in enumerate(Instruction_Set_Conjugated):
            for gate_command in individual_instructions:
                L = len(gate_command)
                if L == 1:
                    gate_command[0] | qubits[idx]
                elif L == 2:
                    C(gate_command[0]) | (qubits[gate_command[1]] , qubits[idx])
        # c-Pi
        Grover_Reflection | qubits

    # R
        
        for idx,individual_instructions in enumerate(Instruction_Set):
            for gate_command in individual_instructions:
                L = len(gate_command)
                if L == 1:
                    gate_command[0] | qubits[idx]
                elif L == 2:
                    C(gate_command[0]) | (qubits[gate_command[1]] , qubits[idx])

    for idx,individual_instructions in enumerate(Instruction_Set_Conjugated):
        for gate_command in individual_instructions:
            L = len(gate_command)
            if L == 1:
                gate_command[0] | qubits[idx]
            elif L == 2:
                C(gate_command[0]) | (qubits[gate_command[1]] , qubits[idx])

    All(Measure) | qubits

    Engine.flush()

    results = 0
    for q in qubits:
        results += int(q)

    if results%2 == 0:
        outcome = 0
    else:
        outcome = 1

    return(
        outcome, M
    )

sig = pi/4
mu = random.uniform(0, 2*pi)
ans_param = 1.23452
rel_val = cos(ans_param)**2
run = 0
Maxruns = 500

Max_M = 10000

while sig > 5*10**-3:
    out, m = circuit_no_ancilla(param = ans_param, sigma = sig, maxM = Max_M)
    mu, sig = ExactUpdate(outcome = out, sigma = sig, mu = mu, M = m)
    # print(abs(np.cos(mu / 2)), sig)
    run += 1
    if run > Maxruns:
        break

print(rel_val)
print(abs(np.cos(mu / 2)))
