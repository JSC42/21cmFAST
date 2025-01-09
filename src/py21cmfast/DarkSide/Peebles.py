LineWidth = 2
FontSize = 18

import numpy as np
from PyLab import Hubble, plt, LCDM_HyRec
import ctypes

def Peebles_Py(z, xe, T):
    Q = 1.602176634E-19
    alpha = 1/137.03599976
    E0 = 13.6 * Q
    kB = 1.38064852E-23
    hb = 1.0545718002693302e-34
    c = 299792458
    me = 9.109383632044565e-31
    np0 = 0.1901567053460595 # number density of H nuclei today

    # Beta2
    x = E0/(kB * T)
    fx1 = np.log(x)
    fx2 = x * np.exp(x/4)
    fx = fx1/fx2
    beta2 = 9.78 * (E0/(2 * np.pi))**1.5 * alpha**2/(me**0.5 * c * hb) * fx

    L2y = 8.227
    xH = np.max([0, 1-xe]) # avoid tinest supplus from xe=1 which can give incorrect Peebles C
    nH = np0 * (1+z)**3 * xH # Technically this is only correct for 21cmFAST xe definition
    H = Hubble(z=z)
    La1 = H * (3 * E0)**3
    La2 = nH * (8*np.pi)**2 * hb**3 * c**3
    La = La1/La2
    # C = (La + L2y)/(La + L2y + beta2) # This can give numerical issue when xe->0
    C = 1 - beta2/(La + L2y + beta2)
    if 1-xe < 1E-50: C= 1 # Analytically we expect C = 1 if xe=1
    if C > 1 or C < 0: 
        C = np.nan
    return C
    
def Peebles_C(z, xe, T, delta = 0):
    '''
    Prepare lib:
    gcc -shared -o Peebles.so Peebles.c; mv Peebles.so test
    '''
    c_lib = ctypes.CDLL('test/Peebles.so')
    Double = ctypes.c_double
    c_function = c_lib.PeeblesFactor
    c_function.argtypes = (Double, Double, Double)
    c_function.restype = Double
    result = c_function(z, xe, T, delta)
    return result

# Testing
z = np.linspace(0, 1600, 1000)
xe, T = LCDM_HyRec(z=z, Use_EoR=1)
c1 = np.zeros(len(z))
c2 = np.zeros(len(z))

for idx, z_ in enumerate(z):
    xe_, T_ = xe[idx], T[idx]
    c1[idx] = Peebles_Py(z=z_, xe = xe_, T=T_)
    c2[idx] = Peebles_C(z=z_, xe = xe_, T=T_)
    
plt.rcParams.update({
    'text.usetex': True,
    'font.family':'Times',
    'text.latex.preamble': r'\usepackage{newtxtext,newtxmath}'})
fig, ax = plt.subplots()
ax.grid(True, which='both', linewidth = 0.3)  # `which='both'` enables major and minor grids
fig.set_size_inches(10, 8)

plt.plot(1+z, c1, 'k', linewidth=LineWidth, label = 'Python')
plt.plot(1+z, c2, '--r', linewidth=LineWidth, label = 'c')
plt.xlabel('$1+z$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$C$',fontsize=FontSize,fontname='Times New Roman')
plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.legend(fontsize=FontSize,loc = 'lower left')
plt.tight_layout()
plt.savefig('/Users/cangtao/Desktop/tmp.pdf')
