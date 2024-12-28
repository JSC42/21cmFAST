LineWidth = 2
FontSize = 18

import matplotlib.pyplot as plt
import h5py
from PyLab import Hubble

Transfer_File = 'data/TransferFunctions.h5'
f = h5py.File(Transfer_File, 'r')
z = f['Axis/z'][:]
f.close()
nz = len(z)

dz = z[0:nz-1] - z[1:nz]
z = z[0:nz-1]
# comoving dist: dx = -c*dz/H
H = Hubble(z=z)
dx = -299792458 * dz/H/3.086E22

plt.rcParams.update({
    'text.usetex': True,
    'font.family':'Times',
    'text.latex.preamble': r'\usepackage{newtxtext,newtxmath}'})

plt.loglog(z, dx, 'k', linewidth=LineWidth, label = 'sin')
plt.xlabel('$z$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('Resolution [cMpc]',fontsize=FontSize,fontname='Times New Roman')
plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
#plt.xlim([11, 60])
plt.tight_layout()
plt.savefig('/Users/cangtao/Desktop/tmp.pdf')
