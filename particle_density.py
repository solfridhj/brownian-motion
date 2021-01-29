import numpy as np
from matplotlib import pyplot as plt

ensemble = np.genfromtxt("many_particles/ensemble.csv", delimiter=',')
time = np.genfromtxt("many_particles/times.csv", delimiter=',')
