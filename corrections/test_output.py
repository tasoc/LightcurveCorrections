import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import glob

globlist = glob.glob('../data/Rasmus/toutput2/*.noisy_detrend')

for sfile in globlist:
    f = np.loadtxt(sfile).T
    plt.plot(f[0], f[1])
    plt.plot(f[0], f[2])
    plt.show()
