#!/usr/bin/env python
import os,sys,glob,time,copy
from math import floor,ceil,log,exp
import numpy as np
import dicom
import matplotlib.pyplot as plt
def main():
  # ==================================================================
  # Configuration
  # ==================================================================
  A = {}
  for l in open('run.conf').readlines():
    v = l.split()  
    try: A[v[0]] = eval(v[1])
    except: A[v[0]] = v[1]
  if len(sys.argv)>1: A.update( eval(sys.argv[1]) )


  print A['DATA_PATH']
  from ROOT import TFile,TH3S,TH3D,TGraph,TVector3,Double

  os.chdir(A['DATA_PATH'])
  try:os.mkdir('out')
  except:pass
  os.chdir('out')
  # ==================================================================
  # Phantom
  # ==================================================================
  voxel_rebin = 2
  # Get CT data
  f = TFile('phantom.root','recreate')
  fCT = glob.glob(A['DATA_PATH']+'/*.dcm')

  fCT.sort()

  CT = dicom.read_file(fCT[-1])
  
  format_str = '%sint%d' % (('u', '')[CT.PixelRepresentation],CT.BitsAllocated)
  N = [0,0,0]
  N[0] = CT.Rows
  N[1] = CT.Columns
  N[2] = len(fCT)
  if not CT.SliceThickness==0: pix = CT.PixelSpacing+[CT.SliceThickness]
  else : 
    try :pix = CT.PixelSpacing + [abs(dicom.read_file(fCT[0]).SliceLocation-dicom.read_file(fCT[1]).SliceLocation)]
    except :pix = CT.PixelSpacing + [abs(dicom.read_file(fCT[0]).ImagePositionPatient[2]-dicom.read_file(fCT[1]).ImagePositionPatient[2])]
    
  
  pos = CT.ImagePositionPatient
  
  # Find the minimum z position
  pos[2] = min([dicom.read_file(dcm).ImagePositionPatient[2] for dcm in fCT] )
  maxZ = max( [dicom.read_file(dcm).ImagePositionPatient[2] for dcm in fCT] )
  minZ = pos[2]
  pix[2] = (maxZ-minZ)/(len(fCT)-1)

  pix=map(float,pix)
  pos=map(float,pos)
  n  = [N[0]/voxel_rebin,N[1]/voxel_rebin,N[2]]
  for i in range(3): pos[i] -= 0.5*pix[i]
                                                                
  h = TH3S('hu','',n[0],pos[0],pos[0]+n[0]*pix[0],n[1],pos[1],pos[1]+n[1]*pix[1],n[2],pos[2],pos[2]+n[2]*pix[2])
  # Cropped phantom
  f.cd()
  print 'Voxel size (mm): original(%.3f,%.3f,%.3f) ' % tuple(pix)
  print 'N voxels: %d x %d x %d' % tuple(N)
  print 'Rebin by a factor : %d' % tuple([voxel_rebin])
  print '           original            '
  for i,I in enumerate(['x','y','z']): print '%s = [ %3.4f , %3.4f ]' % (I,0.1*pos[i],0.1*(pos[i]+N[i]*pix[i]))


  
  a = np.zeros((n[2],n[1],n[0]),'int16')
  for dcm in fCT:
    CT = dicom.read_file(dcm)
    z = h.GetNbinsZ()-h.GetZaxis().FindBin(CT.ImagePositionPatient[2])-1
    if(z==N[2]):z=z-1
    c    = np.fromstring(CT.PixelData,np.dtype(format_str)).reshape((N[1],N[0]))
    c    = rebin(c,(N[1]/voxel_rebin,N[0]/voxel_rebin))
    a[z] = c
      

  a -= 2048
  #a = a*int(dicom.read_file(fCT[0]).RescaleSlope)+int(dicom.read_file(fCT[0]).RescaleIntercept)
  a /= 16  
  a = np.clip(a,-128,128)
  b = np.zeros((n[2]+2,n[1]+2,n[0]+2),'int16')
  b[1:-1,1:-1,1:-1] = a #np.clip(a,1,128+72)
  b = b.flatten()
  h.Set(len(b),b)
  h.Write()
  f.Close()

  
  print A['DATA_PATH']+'/out/phantom.root created: check the phantom'
  sys.exit()

def rebin(a, shape):
  sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
  return a.reshape(sh).mean(-1).mean(1)

if __name__=="__main__":
   main()

