from qutip import *
import matplotlib.pyplot as plt
import numpy as np
import math

# System Parameters
#-----------------------------------
Nc = 2                      # Number of cavity states
Nm = 2                      # Number of mech states
Nq = 2                      # Number of qubit states
wm = 0.05                   # Mechanical mode frequency (GHZ)   (Previously used 1 GHz)
v = 0.25                       # Rate of change of coupling strength
s = 3                       # Start position of rise of coupling strength
g = 20                      # gamma

# Operators
#-----------------------------------
a = tensor(destroy(Nc), qeye(Nm), qeye(Nq))
b = tensor(qeye(Nc), destroy(Nm), qeye(Nq))
q = tensor(qeye(Nc), qeye(Nm), destroy(Nq))

psi1 = tensor(basis(Nc,1),basis(Nm,0),basis(Nq,0))
psi3 = tensor(basis(Nc,0),basis(Nm,0),basis(Nq,1))

psi0 = tensor(basis(Nc,1),basis(Nm,0),basis(Nq,0))      # Initially in mode |a b q> = |1 0 0>


#--------Vitanov Coupling-----------------
def G_1 (t,args):                                  # Coupling Strength
    theta = (math.pi/2)*1/(1+np.exp(-v*(t-s/v)))
    G1 = np.sin(theta)
    return G1

def G_2 (t,args):
    theta = (math.pi/2)*1/(1+np.exp(-v*(t-s/v)))
    G2 = np.cos(theta)
    return G2

def dG_1 (t,args):
    theta = (math.pi/2)*1/(1+np.exp(-v*(t-s/v)))
    dtheta = (math.pi/2)*v*np.exp(-v*(t-s/v))/(1+np.exp(-v*(t-s/v)))**2
    G1 = np.sin(theta)
    dG1 = np.cos(theta)*dtheta
    return dG1

def dG_2 (t,args):
    theta = (math.pi/2)*1/(1+np.exp(-v*(t-s/v)))
    dtheta = (math.pi/2)*v*np.exp(-v*(t-s/v))/(1+np.exp(-v*(t-s/v)))**2
    G2 = np.cos(theta)
    dG2 = -np.sin(theta)*dtheta
    return dG2

def G(t,args):
    G = -g*(G_1(t,0)*dG_2(t,0) - dG_1(t,0)*G_2(t,0))/(G_1(t,0)**2 + (g*G_2(t,0))**2)
    return G


# Hamiltonian
#--------------------------------------
Hom1 = a.dag()*b+a*b.dag()
HQ = g*(b*q.dag()+b.dag()*q)
H1 = 1j*psi1*psi3.dag() - 1j* psi3*psi1.dag()
H = [[H1, G], [Hom1, G_1], [HQ, G_2]]     

tlist = np.linspace(0,50,500)                        # Use range according to v

e_ops = [a.dag()*a, b.dag()*b, q.dag()*q]

result = mesolve(H, psi0, tlist, [], e_ops)  

# Plots
#--------------------------------------
fig, ax = plt.subplots(figsize=(8,7.3))                     #  figsize = (width, height)
ax.plot(tlist, result.expect[0], label="Photon", linewidth=2)
ax.plot(tlist, result.expect[1], label="Phonon")
ax.plot(tlist, result.expect[2], label="Qubit", linewidth=2)
ax.set_ylim([-0.03,1.03])
legend = ax.legend(loc='center right', shadow=True, fontsize = 'x-large') # Legend box position, center/ upper
ax.set_xlabel('Time (ns)', fontsize=28)
ax.set_ylabel('Population', fontsize=28)
ax.set_title('(d)',loc='right', fontsize=28)            # Put figure no. according to v 
ax.tick_params(axis='both', which='major', labelsize=20)
plt.show()