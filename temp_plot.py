import numpy as np
import matplotlib.pyplot as plt

mfs = np.genfromtxt("./single_runs/mf_evos.csv")
cs = np.genfromtxt("./single_runs/conflict_evos.csv")
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(mfs, cs, c=np.arange(cs.shape[0]))
plt.show(fig)
