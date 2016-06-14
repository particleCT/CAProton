from ROOT import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc,rcParams
import sys
#f=TFile("200MeVSlab.root")
f=TFile(sys.argv[1])
phant = TFile(sys.argv[2])

h = phant.Get("hu")

n = [h.GetNbinsZ()+2,h.GetNbinsY()+2,h.GetNbinsX()+2]
Zmin  = h.GetZaxis().GetXmin()
Zmax  = h.GetZaxis().GetXmax()
Zbins = h.GetNbinsZ()
Zaxis = np.linspace(Zmin,Zbins,Zmax)

Xmin  = h.GetXaxis().GetXmin()
Xmax  = h.GetXaxis().GetXmax()

Ymin  = h.GetYaxis().GetXmin()
Ymax  = h.GetYaxis().GetXmax()

n2 = (n[0])*(n[1])*(n[2])
d = h.GetArray()
d.SetSize(n2)
HU = np.array(d,np.int64)
HU = np.reshape(HU,n)

from matplotlib.patches import FancyArrowPatch
phase = f.Get("phase")
for event in phase:
    print phase.x0,phase.y0,phase.z0,phase.Einit,phase.Estop, phase.Einit-phase.Estop
    fig, ax = plt.subplots()
    tracks_Z = phase.tracks_Z
    tracks_X = phase.tracks_X
    tracks_Y = phase.tracks_Y
    ELoss    = phase.ELoss
    print np.sum(ELoss), np.array(ELoss).tolist()
    Zid      = np.abs(Zaxis-phase.z0).argmin()
    print Zid
    plt.plot(tracks_X,tracks_Z,"r-",label="Proton Path")
    plt.imshow(HU[Zid],extent=[Xmin,Xmax,Ymin,Ymax],cmap="gray")
    plt.show()
    

