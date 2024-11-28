LineWidth = 2
FontSize = 18

# Compute deposition efficiencies from DM annihilation, this largely follows 10.1088/1475-7516/2022/03/012
import numpy as np
import h5py
import matplotlib.pyplot as plt
from PyLab import Hubble, TimeNow, os
import py21cmfast as p21c
from joblib import Parallel, delayed
import ctypes

def Load_Data():
    '''
    Load required datasets
    Note that axis doc in TransferFunctions.h5 are written for MatLab
    '''

    # 1. Boost factor
    BoostData = np.load('data/BoostTemplates.npz')
    z = BoostData['z']
    HMF = BoostData['HMF']
    PS = BoostData['PS']
    B = BoostData['B']
    
    # 2. Transfer function
    Transfer_File = 'data/TransferFunctions.h5'
    f = h5py.File(Transfer_File, 'r')
    TF = f['T_Array'][:]
    TF_z = f['Axis/z'][:]
    TF_E = f['Axis/Kinetic_Energy_GeV'][:]
    f.close()
    LnZp = np.log(1+TF_z)
    dLnZp = LnZp[0:len(LnZp)-1] - LnZp[1:len(LnZp)]
    dLnZp = np.sum(dLnZp)/len(dLnZp)
    
    # Organize data
    BoostData = {'z': z, 'B':B, 'HMF' : HMF, 'POWER_SPECTRUM':PS}
    TransferData = {'z': TF_z, 'E': TF_E, 'TF': TF, 'dLnZp': dLnZp}
    r = {'BoostData': BoostData, 'TransferData': TransferData}
    return r

DataSet = Load_Data()

def GetEFF_Kernel(idxE=20, idxz=0, spec = 'e', channel='HIon'):
    '''
    Compute deposition efficiency for annihilating DM
    -- inputs --
    channel:
        HIon
        LyA
        Heat
    '''
    HMF = 0
    PS=2
    UseBoost = 1
    CalibrateB = 0
    # Energy is in GeV
    m = 0.5109989E-3 if spec == 'e' else 0
    Ek = DataSet['TransferData']['E'][idxE]
    z = DataSet['TransferData']['z'][idxz]
    B = 10**np.interp(x=z, xp = DataSet['BoostData']['z'], fp = np.log10(DataSet['BoostData']['B'][HMF, PS, :]), left=np.nan, right=0)
    H = Hubble(z)
    sid = 0 if spec == 'y' else 1
    if channel == 'HIon':
        cid = 0
    elif channel == 'LyA':
        cid = 2
    elif channel == 'Heat':
        cid = 3
    else:
        raise Exception('Wrong choice of channel')
    if not UseBoost: B = 1
    Prefix = H * Ek / (B * (m + Ek) * (1+z)**3)

    # Now the next bit
    zi = DataSet['TransferData']['z']
    
    Bi = 10**np.interp(x=zi, xp = DataSet['BoostData']['z'], fp = np.log10(DataSet['BoostData']['B'][HMF, PS, :]), left=np.nan, right=0)
    if CalibrateB == 1:
        # By logic of 0th order integration, using either the left or right points are both ok
        dx = DataSet['TransferData']['dLnZp']
        x0 = np.log(1+zi)
        x1 = x0+dx
        zi_ = np.exp(x1)
        Bi = 10**np.interp(x=zi_, xp = DataSet['BoostData']['z'], fp = np.log10(DataSet['BoostData']['B'][HMF, PS, :]), left=np.nan, right=0)
    if not UseBoost: Bi = 1
    Hi = Hubble(zi)
    T = DataSet['TransferData']['TF'][:, idxE, idxz, cid, sid]
    fun = Prefix*(1+zi)**3 * Bi * T / Hi
    fc = np.sum(fun)
    w = fun/fc
    return w

x = np.arange(0, len(DataSet['TransferData']['z']))
w = GetEFF_Kernel()
plt.plot(x, w, '+k', linewidth=LineWidth, label = 'sin')
plt.xlabel('index',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('weigh',fontsize=FontSize,fontname='Times New Roman')
plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
#plt.xlim([11, 60])
plt.tight_layout()
plt.savefig('/Users/cangtao/Desktop/tmp.pdf')
