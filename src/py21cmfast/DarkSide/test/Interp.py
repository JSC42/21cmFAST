
# Check that the c interpolation in my 21cmFAST agrees with python table, no need to define a new python interp I think
import numpy as np
from PyLab import plt, Find_Index
import ctypes, tqdm
EFF_Data = np.load('../data/EEE_Table.npz')

def EFF_Interp_c(Ek, z, particle, channel, HMF, PS, UseBoost):
    '''
    C function used in Interp_EFF of 21cmFAST DarkSide.h
    '''
    m = Ek+0.5109989E-3 if particle else Ek
    c_lib = ctypes.CDLL('/Users/cangtao/FileVault/GitHub/21cmFAST/src/py21cmfast/DarkSide/c_products/DarkSide.so')
    Double = ctypes.c_double
    Int = ctypes.c_int
    c_function = c_lib.Interp_EFF
    c_function.argtypes = (Double, Double, Int, Int, Int, Int, Int)
    c_function.restype = Double
    result = c_function(z, m, particle, channel, HMF, PS, UseBoost)
    return result

def CompareTable(eid, particle, channel, HMF, PS, UseBoost):
    fc_ax = EFF_Data['fc']
    Ek = EFF_Data['Ek'][eid]
    z = EFF_Data['z']
    fp = fc_ax[UseBoost, particle, channel, HMF, PS, eid, :]
    fc = []
    for z_ in z:
        fc.append(EFF_Interp_c(Ek, z_, particle, channel, HMF, PS, UseBoost))
    return z, fp, fc

def EFF_Interp_py(Ek, z, particle, channel, HMF, PS, UseBoost):
    '''
    Actually u know what let's write a python version real quick that takes any Ek
    follows the same logic as Interp_EFF
    '''
    fc_ax = EFF_Data['fc'][UseBoost, particle, channel, HMF, PS, :, :] # [idxE, idxz]
    fc_ax = np.log10(fc_ax)
    Ek_ax = EFF_Data['Ek']
    z_ax = EFF_Data['z']
    
    # do z
    id1 = Find_Index(x = z, x_axis = z_ax)
    id2 = id1+1
    x1 = np.log10(1+z_ax[id1])
    x2 = np.log10(1+z_ax[id2])
    f1 = fc_ax[:, id1]
    f2 = fc_ax[:, id2]
    x = np.log10(1+z)
    f = (f2 - f1)*(x-x1)/(x2 - x1) + f1
    
    # do E
    id1 = Find_Index(x = Ek, x_axis = Ek_ax)
    id2 = id1+1
    x1 = np.log10(Ek_ax[id1])
    x2 = np.log10(Ek_ax[id2])
    f1 = f[id1]
    f2 = f[id2]
    x = np.log10(Ek)
    f = (f2 - f1)*(x-x1)/(x2 - x1) + f1
    r = 10**f
    return r

def GetRandomE():
    rd = np.random.random(40)
    Ek_ax = np.log10(EFF_Data['Ek'])
    x1 = np.min(Ek_ax)
    x2 = np.max(Ek_ax)
    dx = x2 - x1
    Ek = x1 + rd*dx
    Ek = 10**Ek
    return Ek

def GetRandomZ():
    rd = np.random.random(63)
    zpax = np.log10(1+EFF_Data['z'])
    x1 = np.min(zpax)
    x2 = np.max(zpax)
    dx = x2 - x1
    zp = x1 + rd*dx
    z = 10**zp - 1
    return z

Ek_Random = GetRandomE()
Z_Random = GetRandomZ()

dif_table = 0
dif_random = 0
count_table = 0
count_random = 0

for eid in tqdm.tqdm(range(40), desc = 'Comparing dif'):
    #for eid in np.arange(0, 40):
    for particle in [0, 1]:
        for channel in [0, 1, 2]:
            for HMF in [0, 1, 2, 3]:
                for PS in [0, 1, 2, 3, 4]:
                    for UseBoost in [0, 1]:
                        # Check Table
                        z, fp, fc = CompareTable(eid, particle, channel, HMF, PS, UseBoost)
                        dif_table += np.sum(np.abs(1-fp/fc))/len(z)
                        count_table +=1
                        # Check random Ek, can 
                        for z_ in Z_Random:
                            Ek = Ek_Random[eid]
                            f1 = EFF_Interp_c(Ek, z_, particle, channel, HMF, PS, UseBoost)
                            f2 = EFF_Interp_py(Ek, z_, particle, channel, HMF, PS, UseBoost)
                            dif_random += np.abs(1-f1/f2)
                            count_random += 1
                        
dif_table = dif_table / count_table
dif_random = dif_random / count_random

print('dif_table=', dif_table)
print('dif_random=', dif_random)

# Plot the final iteration
plt.loglog(z, fc, 'k')
plt.loglog(z, fp, '+r')
plt.show()
