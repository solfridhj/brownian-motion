import numpy as np
from matplotlib import pyplot as plt

time = np.genfromtxt("many_particles/time.csv", delimiter=',')
ensemble = np.genfromtxt("many_particles/ensemble.csv", delimiter=',')

ensemble_r2 = np.genfromtxt("many_particles/ensemble_r2.csv", delimiter=',')

# TO GET A PARTICE, AT ALL TIMES [:, 0]
# TO GET A TIME, FOR ALL PARTICLES [0, :]

time_index = 1

time_val = time[time_index]
print(time_val)

step = 1000


for i in range(0, len(time), step):
    print(i)
    ens = ensemble[i, ]
    ens_2 = ensemble_r2[i, ]
    weig = np.ones_like(ens) / len(ens)
    weig_2 = np.ones_like(ens_2) / len(ens_2)
    # Array of positions for the given time
    time_val = time[i]
    plt.hist(ensemble[i, :], weights=weig, bins=200, label=r"$r_2$")
    plt.hist(ensemble_r2[i, :], weights=weig_2, bins=200, label=r"$r_2$")
    plt.title(r"$t = $" + f"{time_val.round(5)}, $\Delta U=80$ eV")
    plt.xlabel(r"$\mu x$ [m]")
    plt.ylabel(r"$n(x, t)$")
    plt.savefig("density_fig/density" + str(i) +"_flash.pdf", dpi=200)
    plt.show()


