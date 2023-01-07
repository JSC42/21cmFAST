from __future__ import print_function 
import numpy as np
import matplotlib.pyplot as plt
from colossus.cosmology import cosmology
from colossus.halo import profile_nfw
from colossus.halo import concentration

from scipy.integrate import quad

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
  # Get concentration param
  # Tutorial: https://bdiemer.bitbucket.io/colossus/_static/tutorial_halo_concentration.html
  c=concentration.concentration(M=m/h, mdef='vir', z=z)
  p_nfw = profile_nfw.NFWProfile(M = m/h, c = c, z = z, mdef = 'vir')
  rho = p_nfw.density(r/h) /h/h
  return rho
  
M=1E15
r = 10**np.arange(0,4,0.001)
rho = r
id=0
for R in r:
  rho[id] = Profile(M,R,0)
  id = id+1

y = 4*np.pi*r*r*rho
mm= np.trapz(y=y,x=r)
print(mm)