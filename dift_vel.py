
import numpy as np
from matplotlib import pyplot as plt


taus = np.genfromtxt("single_particle/taus_r2.csv", delimiter=',')[:-2]

# erroneously divided by  3000 when computing, so multiply to get rid of (took 40 min. to run)
velocity = np.genfromtxt("single_particle/vel_r2.csv", delimiter=',')[:-2]*3000

# Convert to non-reduced time
deltaU = 80*1.602176565*10**(-19)
L = 20*10**(-6)

eta = 10**(-3)
r_1 = 3*12*10**(-9)

gamma = 6*np.pi*eta*r_1
omega = deltaU/(gamma*L*L)

t_end = omega/3000
print(omega)

# Getting velocity in non-reducedtime
velocity = velocity/t_end
plt.plot(taus, velocity)
plt.xlabel(r"$\tau$ [s]")
plt.ylabel(r"$\langle v \rangle$ [$\mu m/s$]")
plt.show()
#plt.savefig("figures/driftvelocity_r2.pdf", dpi=200)

max_index = np.argmax(velocity)
print(np.max(velocity))
print(max_index)

tau_max = taus[max_index]
print(tau_max)