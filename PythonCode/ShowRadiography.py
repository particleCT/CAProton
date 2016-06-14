import rootplot
import rootplot.root2matplotlib as r2m
from ROOT import *
from matplotlib import pyplot as plt
import sys
import numpy as np
from numpy import *
from scipy.signal import argrelextrema
f=TFile(sys.argv[1])
A=f.Get("pRad1")
histA = r2m.Hist2D(A)
histA = np.array(histA.content)

histA = np.clip(histA,160,300)[120:830,200:800]
#histA =np.rot90(histA,2)
#plt.xlim([100,400])
#plt.ylim([120,450])
histA =np.rot90(histA,2)
plt.imshow(histA,cmap="gray",extent=[-16.5,16.5,-16.5,8.75])

plt.colorbar()
plt.show()
