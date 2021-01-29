
import numpy as np
from matplotlib import pyplot as plt

# PARAMETERS
k_BT = 26
deltaU = 60* k_BT

print("HALLO")
tau = str(0)


def single_particle_plot():
    x = np.genfromtxt("single_particle/x_val.csv", delimiter=',')
    t = np.genfromtxt("single_particle/t_val.csv", delimiter=',')
    x_mean = np.genfromtxt("single_particle/x_mean.csv", delimiter=',')

    # Plotting the potential
    x_vals = np.genfromtxt("potential/x_data.csv", delimiter=',')
    U_vals = np.genfromtxt("potential/potential_data.csv", delimiter=',')

    #pot = np.genfromtxt("potential/potential_data.csv", delimiter=',')
    #x_dat = np.genfromtxt("potential/x_data.csv", delimiter=',')

    #pot = pot[:-1]
    #x_dat = x_dat[:-1]

    #x_1 = np.genfromtxt("single_particle/x_val_1.csv", delimiter=',')
    #t_1 = np.genfromtxt("single_particle/t_val_1.csv", delimiter=',')

    plt.plot(x, t, linewidth=0.8)
    plt.plot(x_vals-1, U_vals*100)
    #plt.plot(x_1, t_1, linewidth=0.8, label=r"$0.1k_BT$")
    plt.plot()
    plt.xlabel(r"$x[\mu m]$")
    plt.ylabel(r"$\hat t$[s]")
    plt.title(r"Trajectory for $\tau = 7.0 s$")
    plt.savefig("figures/trajectory_tau_flash_7.pdf", dpi=200)
    plt.show()
    print("Difference: " + str(x[-2] - x[0]))
    print(x[-2])


single_particle_plot()



