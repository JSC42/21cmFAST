from __future__ import print_function 
import os
import numpy as np
import matplotlib.pyplot as plt
from colossus.cosmology import cosmology
from colossus.halo import profile_nfw
from colossus.halo import concentration
from scipy.integrate import quad
import h5py
import math

def Profile(m,r,z):
  '''Get halo profile
  inputs:
    m: halo mass in msun
    r: radius in kpc
    z: redshift
  return value:
    rho_dm in msun/kpc^3
  '''
  cosmo = cosmology.setCosmology('planck18')
  h = 0.6766
  Mh=m/h
  rh=r/h
  # Get concentration param
  # Tutorial: https://bdiemer.bitbucket.io/colossus/_static/tutorial_halo_concentration.html
  c=concentration.concentration(M = Mh, mdef='vir', z=z)
  p_nfw = profile_nfw.NFWProfile(M = Mh, c = c, z = z, mdef = 'vir')
  rho = p_nfw.density(rh)*h**2
  return rho

def Solve(f,lx1,lx3,dlx):
  Go=1
  Lx1=lx1
  Lx3=lx3
  Count=1
  while Go:
    Lx2=(Lx1+Lx3)/2
    x1=10**Lx1
    x2=10**Lx2
    x3=10**Lx3
    f1=f(x1)
    f2=f(x2)
    f3=f(x3)
    s1=np.int32(np.sign(f1))
    s2=np.int32(np.sign(f2))
    s3=np.int32(np.sign(f3))
    if s1==s3:
      # print('Crash immenent, debug info: f1=',f1,' f3=',f3)
      raise Exception('Solution not found')
    if s1 == s2:
      Lx1 = Lx2
    elif s2 == s3:
      Lx3 = Lx2
    if Lx3-Lx1 < dlx:
      Go=0
      return 10**Lx1
    Count = Count+1
    if Count > 100:
      raise Exception('Solution not found and took too long')

def FindMass(m,R,z):
  r = 10**np.arange(-10.2,math.log10(R),0.01)
  rho=Profile(m,r,z)
  y=4*np.pi*r**2*rho
  M = np.trapz(y=y,x=r)
  return M

def GetRmax(m,z,method):
  '''Get halo profile
  inputs:
    m: halo mass in msun
    z: redshift
  return value:
    rho_dm in msun/kpc^3
  '''
  cosmo = cosmology.setCosmology('planck18')
  h = 0.6766
  Mh=m/h
  c=concentration.concentration(M = Mh, mdef='vir', z=z)
  p_nfw = profile_nfw.NFWProfile(M = Mh, c = c, z = z, mdef = 'vir')
  if method == 1:
    fun = lambda x: p_nfw.enclosedMassInner(x) - m
  else:
    fun = lambda x: FindMass(m,x,z) - m
  Rmax=Solve(fun,-10,25,0.01)
  return Rmax

def GetConvergedProfile(m,z):
  """
  m: in msun
  return value:
  r and rho in SI unit
  """
  Rmax=GetRmax(m,z,2)
  r = 10**np.arange(-10,math.log10(Rmax),0.001)
  rho=Profile(m,r,z)
  kpc = 3.086E19
  msun=1.988E30
  r=r*kpc
  rho=rho*msun/kpc**3
  return r,rho

def IntegrateProfile(m,z):
  h=0.6766
  OmC=0.261
  RhoCr=1.879E-26 *h**2
  print(m,z)
  #try:
  r,rho = GetConvergedProfile(m,z)
  fun=r**2 * rho**2
  return np.trapz(y=fun,x=r)
  #except:
  #  print("Sth went wrong!")
  #  return 0


# ---- Start the run ----
M1=1E5
z1=25.25
a=IntegrateProfile(M1,z1)
print(a)

