from qutip import *
import matplotlib.pyplot as plt
import numpy as np
import math

# System Parameters
#-----------------------------------
Nc = 100                      # Number of cavity states
Nm = 2                      # Number of mech states
Nq = 2                      # Number of qubit states
wm = 1                      # Mechanical mode frequency (GHZ)
v = 0.5                       # Rate of change of coupling strength
s = 3                       # Start position of rise of coupling strength
g = 20                      # gamma


# Operators
#-----------------------------------
a = tensor(destroy(Nc), qeye(Nm), qeye(Nq))
b = tensor(qeye(Nc), destroy(Nm), qeye(Nq))
q = tensor(qeye(Nc), qeye(Nm), destroy(Nq))

psi0 = tensor(basis(Nc,3),basis(Nm,0),basis(Nq,0))      # Initially in mode |a b q> = |1 0 0>


#--------Vitanov Coupling-----------------
def G_1 (t, args):                                  # Coupling Strength
    theta = (math.pi/2)*1/(1+np.exp(-v*(t-s/v)))
    G1 = np.sin(theta)
    return G1

def G_2 (t,args):
    theta = (math.pi/2)*1/(1+np.exp(-v*(t-s/v)))
    G2 = np.cos(theta)
    return G2


# Hamiltonian
#--------------------------------------
H0 = wm*(b.dag()*b+0.5*q.dag()*q)
Hom1 = a.dag()*b+a*b.dag()
HQ = g*(b*q.dag()+b.dag()*q)
H = [[Hom1, G_1], [HQ, G_2], H0]

tlist = np.linspace(0,50,300)                       

e_ops = [a.dag()*a, b.dag()*b, q.dag()*q]

result = mesolve(H, psi0, tlist, [], e_ops)  

# Plots
#--------------------------------------
fig, ax = plt.subplots(figsize=(8,7.2))                     #  figsize = (width, height)
ax.plot(tlist, result.expect[0], label="Photon", linewidth=2.0)
ax.plot(tlist, result.expect[1], label="Phonon")
ax.plot(tlist, result.expect[2], label="Qubit", linewidth=2.0)
# ax.set_ylim([-0.03,1.03])
legend = ax.legend(loc='upper right', shadow=True, fontsize = 'xx-large') # Legend box position, center/ upper
ax.set_xlabel('Time (ns)', fontsize=28)
ax.set_ylabel('Population', fontsize=28)
ax.set_title('(d)',loc='right', fontsize=28)                       # Put figure no. according to v 
ax.tick_params(axis='both', which='major', labelsize=20)
plt.show()