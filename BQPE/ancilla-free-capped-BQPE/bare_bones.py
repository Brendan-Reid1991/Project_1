import numpy as np
from numpy import pi, exp, cos, sin

def Update(M, outcome, mu, sigma):
    Expectation = mu - (
        (1-2*outcome)*M*sigma**2 * sin(M*mu)
        ) / (
            exp(M**2 * sigma**2 / 2) + (1-2*outcome)*cos(M*mu)
            )

    VarNum = exp(M**2 * sigma**2) + (
        0.5*(2*outcome - 1) * (
            2*exp(M**2 * sigma**2 / 2)*(M**2 * sigma**2 - 2)*cos(M*mu) + (2*outcome - 1)*(
                1 - (2 * M**2 * sigma**2) + cos(2*M*mu)
            )
        )
        )

    VarDenom = (
        exp(M**2 * sigma**2 / 2) + (1 - 2*outcome)*cos(M*mu)
        )**2
    
    Variance = sigma**2 * (VarNum / VarDenom)
    
    Std = np.sqrt(Variance)
    return(Expectation, Std)

##########################################################################################
##########################################################################################
##########################################################################################

def Ancilla_Based(phi = 0, mu = 0, sigma = 0):

    runs = 0
    Max_Runs = 10**4

    while sigma > 0.005:

        M = max(1, np.round(1/sigma))
        P0 = 1/2*(1+cos(M * phi))

        x = np.random.uniform(0, 1)
        if x < P0:
            out = 0
        else:
            out = 1
        
        mu, sigma = Update(M = M, outcome = out, mu = mu, sigma = sigma)
         
        
        runs += 1
        if runs > Max_Runs:
            break
    if runs > Max_Runs:
        return(1000,M, runs)
    else:
        return(mu,M, runs)

##########################################################################################
##########################################################################################
##########################################################################################

def Ancilla_Free(phi = 0, mu = 0, sigma = 0):

    runs = 0
    Max_Runs = 10**4

    while sigma > 0.005:

        M = max(1, round(0.5/sigma))
        if M%2 == 0 and M > 1:
            pass
        else:
            M = max(1, M - 1)
        P0 = 1/2*(1 + cos(M * phi))

        x = np.random.uniform(0, 1)
        if x < P0:
            out = 0
        else:
            out = 1
        
        mu, sigma = Update(M = M, outcome = out, mu = mu, sigma = sigma)

        runs += 1
        if runs > Max_Runs:
            break
    if runs > Max_Runs:
        return(1000, M, runs)
    else:
        return(mu, M, runs)

##########################################################################################
##########################################################################################
##########################################################################################

from progress.bar import ShadyBar

Phi_Range = np.linspace(-pi, pi, 100)

average_over = 100

bar = ShadyBar('Processing...', max = len(Phi_Range)*average_over, suffix = "%(percent).2f%%")

Ab = []  ## P0 = 1/2*[1+Cos(M * phi)]
Af = []   ## P0 = 1/2*[1+Cos(2*M * phi)]

Sanscilla_M = []
Concilla_M = []

Sanscilla_Runs = []
Concilla_Runs = []

for phi in Phi_Range:
    temp = [[], []]
    Ms = [[], []]
    Runs = [[], []]
    counter = 0
    while counter < average_over:
        mu_Af, final_M_S, runs_S = Ancilla_Free(phi = phi, mu = np.random.uniform(-pi, pi), sigma = pi/4)

        mu_Ab, final_M_C, runs_C = Ancilla_Based(phi = phi, mu = np.random.uniform(-pi, pi), sigma = pi/4)

        if mu_Af == 1000 or mu_Ab == 1000:
            pass
        else:
            temp[0].append(abs(cos(mu_Ab/2)))
            temp[1].append(abs(cos(mu_Af/2)))

            Ms[0].append(final_M_S)
            Ms[1].append(final_M_C)

            Runs[0].append(runs_S)
            Runs[1].append(runs_C)

            bar.next()  
            counter += 1

    Ab.append(
        np.median(temp[0])
    )
    Af.append(
        np.median(temp[1])
    )

    Sanscilla_Runs.append(np.mean(Runs[0]))
    Concilla_Runs.append(np.mean(Runs[1]))

    Sanscilla_M.append(np.mean(Ms[0]))
    Concilla_M.append(np.mean(Ms[1]))

bar.finish()

np.savetxt("Ancilla_Based_AgainstPhi.txt", np.asarray(Ab))
np.savetxt("Ancilla_Free_AgainstPhi.txt", np.asarray(Af))

np.savetxt("Ancilla_Based_AgainstPhi_M.txt", np.asarray(Concilla_M))
np.savetxt("Ancilla_Free_AgainstPhi_M.txt", np.asarray(Sanscilla_M))

np.savetxt("Ancilla_Based_AgainstPhi_Runs.txt", np.asarray(Concilla_Runs))
np.savetxt("Ancilla_Free_AgainstPhi_Runs.txt", np.asarray(Sanscilla_Runs))

import matplotlib.pyplot as plt
plt.plot(Phi_Range/pi, cos(np.asarray(Phi_Range)/2), linewidth = 1.5, color = "grey", label = "True")
plt.plot(Phi_Range/pi, Ab, linewidth = 1.5, label = 'Concilla', linestyle = 'dotted')
plt.plot(Phi_Range/pi, Af, linewidth = 1.5, label = 'Sanscilla', linestyle = 'dashed')
plt.xlabel(r"$\phi/\pi$")
plt.ylabel(r"$\cos(\mu/2)$")
plt.legend()
plt.grid(True)
plt.savefig('Ab_vs_Af.pdf')
# plt.show()

plt.clf()

plt.plot(Phi_Range/pi, cos(np.asarray(Phi_Range)/2), linewidth = 1.5, color = "grey", label = "True")
plt.plot(Phi_Range/pi, Af, linewidth = 2.0, label = 'Sanscilla')
plt.xlabel(r"$\phi/\pi$")
plt.ylabel(r"$\cos(\mu/2)$")
plt.legend()
plt.title("Ancilla free")
plt.grid(True)
plt.savefig('Af.pdf')
# plt.show()

plt.clf()

plt.plot(Phi_Range/pi, cos(np.asarray(Phi_Range)/2), linewidth = 1.5, color = "grey", label = "True")
plt.plot(Phi_Range/pi, Ab, linewidth = 2.0, label = 'Concilla')
plt.xlabel(r"$\phi/\pi$")
plt.ylabel(r"$\cos(\mu/2)$")
plt.legend()
plt.title("Ancilla based")
plt.grid(True)
plt.savefig('Ab.pdf')
# plt.show()

plt.clf()

plt.plot(Phi_Range/pi, Concilla_M, linewidth = 1.5, label = 'Concilla')
plt.plot(Phi_Range/pi, Sanscilla_M, linewidth = 1.5, label = 'Sanscilla')
plt.xlabel(r"$\phi/\pi$")
plt.ylabel(r"$M$")
plt.legend()
plt.grid(True)
plt.savefig('M_Sans_Vs_Con.pdf')

plt.clf()

plt.plot(Phi_Range/pi, Concilla_Runs, linewidth = 1.5, label = 'Concilla')
plt.plot(Phi_Range/pi, Sanscilla_Runs, linewidth = 1.5, label = 'Sanscilla')
plt.xlabel(r"$\phi/\pi$")
plt.ylabel(r"Runs")
plt.legend()
plt.grid(True)
plt.savefig('Runs_Sans_Vs_Con.pdf')