from ROOT import *
import numpy as np
import matplotlib.pyplot as plt
import sys

phantom = TFile(sys.argv[1],"update")

h = phantom.Get('hu')
n = [h.GetNbinsX(),h.GetNbinsY(),h.GetNbinsZ()]
Nn = [n[2]+2,n[1]+2,n[1]+2]
n2 = (n[0]+2)*(n[1]+2)*(n[2]+2)
d = h.GetArray()
d.SetSize(n2)
pf = np.array(d,np.float64)
pf = pf.reshape(Nn)

f = TFile(sys.argv[2],'update')
h = f.Get('Edep')
d = h.GetArray()
d.SetSize(n2)
D = np.array(d,np.float64)
D = D.reshape(Nn)


D   = np.rot90(D,1).T
pf  = np.rot90(pf,1).T


for x in range(250,300):
    PhantomCT = pf[x][::-1,:]
    Dose      = D[x][::-1,:]
    plt.imshow(PhantomCT,cmap='gray',extent=[-150,143,-135,85.0])
    plt.imshow(Dose,cmap='hot',alpha=0.4,extent=[-150,143,-135,85.0])
    plt.show()

