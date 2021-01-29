import numpy as np
from matplotlib import pyplot as plt
from scipy import integrate

# PARAMETERS

# Boltzman distribution
kb_T = (26*10**(-3))
#kb_T = (0.26)
U = np.linspace(0, 1)
U_kbt = U/kb_T
deltaU_kbt = 10

def boltz(U):
    # return np.exp(-U/k_BT)/(k_BT*1.602176565*10**(-19)*(1-np.exp(-deltaU/k_BT)))
    return np.exp(-U_kbt) / (1 - np.exp(-deltaU_kbt))*1/kb_T

p = boltz(U)

U_avg_01 = np.genfromtxt("boltzmann/U_avg.csv", delimiter=',')
U_10 = np.genfromtxt("boltzmann/U_avg_10.csv", delimiter=',')

# To get normalized

weights_01 = np.ones_like(U_avg_01) / len(U_avg_01)
weights_10 = np.ones_like(U_10) / len(U_10)
plt.hist(U_10, weights=weights_10, bins=1000, label=r"$\Delta U = 10 k_BT$")
plt.hist(U_avg_01, weights=weights_01, bins=1000, label=r"$\Delta U = 0.1k_BT$")

plt.ylabel(r'$p(\hat U)$')
plt.xlabel(r"$\hat U_r$")
plt.title("Distribution of occupied potential energy")
plt.legend()
plt.savefig("figures/boltzmann_dist.pdf", dpi=200)
#plt.show()

