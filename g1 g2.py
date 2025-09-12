import matplotlib.pyplot as plt
import numpy as np
import math

v = 0.25                   # Rate of change of coupling strength
s = 3                      # Start position of rise of coupling strength

#--------Vitanov Coupling-----------------
def G_1v (t):                                  # Coupling Strength
    theta = (math.pi/2)*1/(1+np.exp(-v*(t-s/v)))
    G1 = np.sin(theta)
    return G1

def G_2v (t):
    theta = (math.pi/2)*1/(1+np.exp(-v*(t-s/v)))
    G2 = np.cos(theta)
    return G2

tlist = np.linspace(0,50,300)                       

# Plots
#--------------------------------------
fig, ax = plt.subplots(figsize=(8,7.3))               #  figsize = (width, height)
ax.plot(tlist, G_1v(tlist).astype(np.float), label="$g_1$", linewidth=2.0)
ax.plot(tlist, G_2v(tlist).astype(np.float), label="$g_2$", linewidth=2.0)
ax.set_ylim([-0.03,1.03])
legend = ax.legend(loc='center right', shadow=True, fontsize = 'xx-large')      # Legend box position, center/ upper
ax.set_xlabel('Time', fontsize=28)
ax.set_ylabel('Coupling Strength', fontsize=28)
ax.set_title('(b)',loc='right', fontsize=28)            # Put figure no. according to v 
ax.tick_params(axis='both', which='major', labelsize=22)
plt.show()