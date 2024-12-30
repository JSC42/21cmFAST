spec=1
dep_channel=1
idxe = 30
LineWidth = 2
FontSize = 18

import ctypes
from PyLab import plt, np
# Check that my python results agrees with HyRec
# In HyRec fc is defined as the fraction of injected kinetic energy so some dif at low E is OK
def EFF_Interp_HyRec(z, Ek, dep_channel, particle):
    '''
    Deposition interpolation in HyRec, adapted from Interpolate_EFF.h
    '''
    dep_channels = [1, 3, 4]
    particles = [1, 2]
    m = 0.0 if particle == 0 else 0.5109989E-3
    c_lib = ctypes.CDLL('/Users/cangtao/Desktop/tmp/HyRec/Interpolate_EFF.so')
    Double = ctypes.c_double
    Int = ctypes.c_int
    c_function = c_lib.Interp_EFF_DM_Annihilation
    c_function.argtypes = (Double, Double, Int, Int)
    c_function.restype = Double
    result = c_function(Ek+m, z, dep_channels[dep_channel], particles[particle])
    return result

def Test_results(idxe, spec, dep_channel):
    '''
    Check that things are working properly, for given E-index, spec, dep_channel, the code gives fc from python (f1) and HyRec (f0) respectively
    '''
    UseBoost=0
    HMF=0
    PS=2
    EffData = np.load('/Users/cangtao/FileVault/GitHub/21cmFAST/src/py21cmfast/DarkSide/data/EEE_Table.npz')
    #EffData = np.load('/Users/cangtao/FileVault/GitHub/21cmFAST/src/py21cmfast/DarkSide/data/EEE_Table_tmp_wrong.npz')
    Ek = EffData['Ek'][idxe]
    z = EffData['z']
    f1 = EffData['fc'][UseBoost, spec, dep_channel, HMF, PS, idxe, :]
    nz = len(z)
    f0 = np.zeros(nz)
    for idx, z_ in enumerate(z):
        f0[idx] = EFF_Interp_HyRec(z=z_, Ek=Ek, dep_channel=dep_channel, particle=spec)
    #print('Ek_GeV/me =', Ek/0.5109989E-3)
    return z, f0, f1

def Compare_dif(idxe, spec, dep_channel):
    z, f0, f1 = Test_results(idxe=idxe, spec=spec, dep_channel=dep_channel)
    f0[f0<1E-100] = 1E-100
    f1[f1<1E-100] = 1E-100
    dif = np.sum(np.abs(1-f0/f1))/len(f1)
    return dif

z, f0, f1 = Test_results(idxe=idxe, spec=spec, dep_channel=dep_channel)

plt.rcParams.update({
    'text.usetex': True,
    'font.family':'Times',
    'text.latex.preamble': r'\usepackage{newtxtext,newtxmath}'})
fig, ax = plt.subplots()

plt.loglog(1+z, f0, 'k', linewidth=LineWidth, label = 'HyRec')
plt.loglog(1+z, f1, '+r', linewidth=LineWidth, label = 'Python')
plt.xlabel('$1+z$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$f_{\mathrm{c}}$',fontsize=FontSize,fontname='Times New Roman')
plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.legend(fontsize=FontSize,loc = 'upper left')
plt.tight_layout()
f0[f0<1E-100] = 1E-100
f1[f1<1E-100] = 1E-100
dif = np.sum(np.abs(1-f0/f1))/len(f1)
print('dif =', dif)
plt.savefig('/Users/cangtao/Desktop/tmp.pdf')

# Finally test the dif between all channels
dif = 0
count = 0
for idxe in np.arange(0, 40):
    for spec in [0, 1]:
        for dep in [0, 1, 2]:
            dif += Compare_dif(idxe, spec, dep)
            count += 1
dif_avg = dif/count
print('Averaged dif between all channels, energy and species:', "{:.3f}".format(dif_avg))
if dif_avg > 1E-2:
    raise Exception('Dif too large, maybe we are using results in which fc is defined as ratio of TOTAL energy (with mass)? '+
                    'HyRec uses fraction of kinetic energy!')
