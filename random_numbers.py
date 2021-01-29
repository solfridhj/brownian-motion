import numpy as np
from matplotlib import pyplot as plt


random_values = np.genfromtxt("random/rand_vals.csv", delimiter=',')

print(random_values)

plt.hist(random_values, bins=1000)
plt.ylabel("Frequency")
plt.xlabel("Random Number")
plt.title(r"$N=1000000$ samples")
#plt.show()
plt.savefig("normal_dist.pdf", dpi=200)
print(np.mean(random_values))
print(np.std(random_values))
