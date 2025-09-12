from qutip import *
import matplotlib.pyplot as plt
import numpy as np
import math

# System Parameters
#-----------------------------------
Nc = 2                      # Number of cavity states
Nm = 2                      # Number of mech states
Nq = 2                      # Number of qubit states
k = 0.005                   # Kappa (Cavity 1 Damping Rate) (GHz)
ym = 5*10**(-5)             # Mechanical Damping Rate (GHz)
kq = 0.005                  # Qubit emission rate (GHz)
n_th = 50                   # Thermal Phonon Number
wm = 0.05                   # Mechanical mode frequency (GHz) (Previously used 1 GHz)
v = 0.75                    # Rate of change of coupling strength
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

# Collapse Operator
#-------------------------------------
cc = np.sqrt(k)*a
cm = np.sqrt(ym*(1.0 + n_th))*b
cp = np.sqrt(ym*n_th)*b.dag()
cq = np.sqrt(kq)*q

c_ops = [cc,cm,cp,cq]

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

result = mesolve(H, psi0, tlist, c_ops, [])  

# Fidelity
#--------------------------------------
F = []
rho_m = []              
S01 = tensor(basis(Nc,0),basis(Nc,1))                # State |a q> = |0 1>
for i in range(len(result.states)):
    rho_m.append((result.states[i]).ptrace([0,2]))   # Taking partial trace keeping the 1st and 3rd components
    F.append(abs(S01.dag()*rho_m[i]*S01))            # F = <01|tr_m(rho)|01>       
  

# Plots
#--------------------------------------
fig, axes = plt.subplots(figsize=(8,7.2))            #  figsize = (width, height)
f_temp = np.real(np.array(list(map(lambda var: var.todense()[0][0], F))))
axes.plot(tlist, list(map(lambda var: var[0][0], f_temp)), linewidth=2)
axes.set_xlabel(r'Time (ns)', fontsize=28)
axes.set_ylabel(r"Fidelity", fontsize=28)
axes.set_title('(b)',loc='right', fontsize=28)                  # Put figure no. according to v 
axes.tick_params(axis='both', which='major', labelsize=20)
plt.ylim(0,1)
plt.show()