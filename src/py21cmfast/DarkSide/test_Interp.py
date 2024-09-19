from ComputeEFF import *

# Check that the c interpolation agrees with python
def EFF_Interp_c(Ek, z, particle, channel, HMF, PS, UseBoost):
    m = Ek+0.5109989E-3 if particle else Ek
    c_lib = ctypes.CDLL('/Users/cangtao/FileVault/GitHub/21cmFAST/src/py21cmfast/DarkSide/c_products/Interp_EFF.so')
    Double = ctypes.c_double
    Int = ctypes.c_int
    c_function = c_lib.Interp_EFF
    c_function.argtypes = (Double, Double, Int, Int, Int, Int, Int)
    c_function.restype = Double
    result = c_function(z, m, particle, channel, HMF, PS, UseBoost)
    return result

def EFF_Interp_all(eid, particle, channel, HMF, PS, UseBoost):
    EFF_Data = np.load('data/EEE_Table.npz')
    fc_ax = EFF_Data['fc']
    Ek = EFF_Data['Ek'][eid]
    z = EFF_Data['z']
    fp = fc_ax[UseBoost, particle, channel, HMF, PS, eid, :]
    fc = []
    for z_ in z:
        fc.append(EFF_Interp_c(Ek, z_, particle, channel, HMF, PS, UseBoost))
    return z, fp, fc

'''
eid = 0
particle = 0
channel =0 
HMF = 0
PS = 2
UseBoost = 0
z, fp, fc = EFF_Interp_all(eid, particle, channel, HMF, PS, UseBoost)
dif = np.sum(np.abs(1-fp/fc))/len(z)
print(dif)
plt.loglog(z, fc, 'k')
plt.loglog(z, fp, '+r')
plt.show()
'''
dif = 0
count = 0
for eid in np.arange(0, 40):
    for particle in [0, 1]:
        for channel in [0, 1, 2]:
            for HMF in [0, 1, 2, 3]:
                for PS in [0, 1, 2, 3, 4]:
                    for UseBoost in [0, 1]:
                        z, fp, fc = EFF_Interp_all(eid, particle, channel, HMF, PS, UseBoost)
                        dif += np.sum(np.abs(1-fp/fc))/len(z)
                        count +=1
                        print('status:', count/(40*2*3*4*5*2))
dif = dif / count
print('dif=', dif)
plt.loglog(z, fc, 'k')
plt.loglog(z, fp, '+r')
plt.show()
