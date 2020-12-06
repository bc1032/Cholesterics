import numpy as np
import pandas as pd
import pylab as pl
import sys
import os
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
lz,lx = 100,100

# data = np.loadtxt('results/eigenvalues9lz100bylx100.dat')
for i in range(0,50):
    data = np.loadtxt("qmod/modulusq%dlz100bylx100.csv" % (i))

    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    print(np.shape(data))
    data = np.reshape(data,(lz,lx))
    #print((data[0,:]))
    #plt.hist2d(data[:,0],data[:,1])#,range=[[0,lz],[0,lx]])#,bins = lz*lx)
    plt.imshow(data[:,],cmap='viridis', vmin = 0.0, vmax = 1)#, vmin=0.0, vmax=0.49)#,  interpolation='nearest')
    #plt.colorbar()
    # plt.savefig("eigen0lz100bylx100second.pdf")
    plt.savefig("vids/mod%dlz100bylx100second.png" % i)

#plt.show()
