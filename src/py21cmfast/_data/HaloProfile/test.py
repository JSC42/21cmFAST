from __future__ import print_function 
import numpy as np
import matplotlib.pyplot as plt
from colossus.cosmology import cosmology
from colossus.halo import profile_nfw

cosmo = cosmology.setCosmology('planck18')

Mvir = 1E15
cvir = 5.0
z = 0.0
p_nfw = profile_nfw.NFWProfile(M = Mvir, c = cvir, z = z, mdef = 'vir')

r = 10**np.arange(0,4,0.01)
rho_m = cosmo.rho_m(z)
rho_nfw = p_nfw.density(r)

plt.figure()
plt.loglog()
plt.xlabel('r(kpc/h)')
plt.ylabel('density / mean')
plt.plot(r, rho_nfw / rho_m, '-', label = 'NFW');
plt.ylim(1E0, 1E7)
plt.legend();

h=0.6766
y = 4 * np.pi * r * r * rho_nfw/h**2
mm= np.trapz(y=y,x=r)
print(mm/Mvir)

