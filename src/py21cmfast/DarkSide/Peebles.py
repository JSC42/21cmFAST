T = 1E4
xe = 0.002
z = 214

import numpy as np
from PyLab import Hubble, plt

# \alpha^{(2)} = \frac{9.78\hbar^2}{m^2_{\rm{e}}c}\left(\frac{\epsilon_0}{k_{\rm{B}}T}\right)^{1/2}\ln\left(\frac{\epsilon_0}{k_{\rm{B}}T}\right)
# \beta = \alpha^{(2)}\left(\frac{m_{\rm{e}} k_{\rm{B}}T}{2 \pi \hbar^2}\right)^{3/2}{\rm{e}}^{-\epsilon_0/k_{\rm{B}}T}
# \Lambda_\alpha=\frac{H(3 \epsilon_0)^3}{n_{\rm{H}}(8 \pi)^2\hbar^3c^3}
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
nH = np0 * (1+z)**3 * (1-xe)
H = Hubble(z=z)
La1 = H * (3 * E0)**3
La2 = nH * (8*np.pi)**2 * hb**3 * c**3
La = La1/La2

C1 = (La + L2y)/(La + L2y + beta2)
# test c
import ctypes
def PeeblesFactor(z_, xe_, Tk_):
    c_lib = ctypes.CDLL('/Users/cangtao/FileVault/GitHub/21cmFAST/src/py21cmfast/DarkSide/Peebles.so')
    Double = ctypes.c_double
    c_function = c_lib.PeeblesFactor
    c_function.argtypes = (Double, Double, Double)
    c_function.restype = Double
    result = c_function(z_, xe_, Tk_)
    return result

C2 = PeeblesFactor(z, xe, T)
dif = (1-C1/C2)
print('C1 =', C1, 'C2 =', C2, 'dif =',dif)
