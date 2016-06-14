### Download the data, if not already on disk #################################

### Read the data in a numpy 3D array #########################################
import mayavi.mlab as mlab
import numpy as np
from ROOT import *
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
from tvtk.util.ctf import ColorTransferFunction
import rootplot.root2matplotlib as r2m
from ROOT import *
from matplotlib import pyplot as plt
from tvtk.util.ctf import PiecewiseFunction
f2=TFile("dose.root")
A=f2.Get("pRad1")
histA = r2m.Hist2D(A)
histA = np.array(histA.content)
histA = np.clip(histA,290,330)
histA = histA[150:850,150:850]
f=TFile("CIRSHead.root")
HU = f.Get("hu")
print HU.GetNbinsX(),HU.GetNbinsY(),HU.GetNbinsZ()
n = [HU.GetNbinsX(),HU.GetNbinsY(),HU.GetNbinsZ()]
n2 = (n[0]+2)*(n[1]+2)*(n[2]+2)
d = HU.GetArray()
d.SetSize(n2)
ctf = ColorTransferFunction()
ctf.add_rgb_point(0.5, 1, 1,1)  # r, g, and b are float
                                   # between 0 and 1

otf = PiecewiseFunction()
otf.add_point(-110, 0.0)
otf.add_point(-105, 1.0)
otf.add_point(128, 1.0)

#D = np.array(d).reshape(282,258,258)
D = np.array(d).reshape(280+2,512+2,512+2)

#D = D[80:270,37:210,64:193]
X, Y ,Z = np.where(D)
mlab.figure(bgcolor=(0, 0, 0), size=(400, 400))
vol = mlab.pipeline.volume(mlab.pipeline.scalar_field(D,opacity=1.0),vmin=-60, vmax=128)
vol.actors.position = [0, -100, 0]
vol._volume_property.set_color(ctf)
vol._volume_property.set_scalar_opacity(otf)
vol._ctf = ctf
vol._otf = otf
vol.update_ctf = True
obj = mlab.imshow(histA,colormap="bone")
obj.actor.orientation = [90, 0, 0] # the required orientation 
obj.actor.position = [95, 194, 65] # the required  position 
obj.actor.scale = [190./700.,153./700. , 0]
obj2 = mlab.imshow(histA,colormap="bone")
obj2.actor.orientation = [0, 0, 0] # the required orientation 
obj2.actor.position = [95, 95, -34] # the required  position 
obj2.actor.scale = [190./700.,153./700. , 0]
#D = ndimage.zoom(D,0.5) 

# Display the data ############################################################


#D = D[np.where(D>-10)]
#mlab.points3d(X,Y,Z,D,colormap="Greys")

#src.spacing = [1, 1, 1.5]
#src.update_image_data = True

mlab.show()

