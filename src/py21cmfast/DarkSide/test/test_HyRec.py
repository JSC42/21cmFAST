from ComputeEFF import *

# Check that my python results agrees with HyRec
def EFF_Interp_HyRec(z, Ek, dep_channel, particle):
    '''
    Deposition interpolation in HyRec, adapted from Interpolate_EFF.h
    '''
    dep_channels = [1, 3, 4]
    particles = [1, 2]
    m = 0.0 if particle == 0 else 0.5109989E-3
    c_lib = ctypes.CDLL('/Users/cangtao/FileVault/GitHub/21cmFAST/src/py21cmfast/DarkSide/test/main.so')
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
    EffData = np.load('data/EEE_Table.npz')
    Ek = EffData['Ek'][idxe]
    z = EffData['z']
    f1 = EffData['fc'][UseBoost, spec, dep_channel, HMF, PS, idxe, :]
    nz = len(z)
    f0 = np.zeros(nz)
    for idx, z_ in enumerate(z):
        f0[idx] = EFF_Interp_HyRec(z=z_, Ek=Ek, dep_channel=dep_channel, particle=spec)
    print('Ek_GeV/me =', Ek/0.5109989E-3)
    return z, f0, f1
'''
spec=0
dep_channel=2
idxe = 39
z, f0, f1 = Test_results(idxe=idxe, spec=spec, dep_channel=dep_channel)

plt.loglog(z, f0, 'k')
plt.loglog(z, f1, '+r')

f0[f0<1E-100] = 1E-100
f1[f1<1E-100] = 1E-100
dif = np.sum(np.abs(1-f0/f1))/len(f1)
print('dif =', dif)
plt.show()
'''