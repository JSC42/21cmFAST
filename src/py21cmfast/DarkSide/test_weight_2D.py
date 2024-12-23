Species = 0
LineWidth = 1.2
FontSize = 9

# What is contribution to fc from different z timesteps? 2D version
import numpy as np
import h5py, time
import matplotlib.pyplot as plt
from PyLab import Hubble

def Load_Data():
    #Load required datasets

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

def Compute_Weigh(idxE=20, idxz=0, spec = 'e', channel='HIon'):
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
    w = fun[idxz]/fc
    return w

zax = DataSet['TransferData']['z']
Eax = DataSet['TransferData']['E']
def Find_idx(z):
    dif = np.abs(z-zax)
    r = np.argmin(dif)
    return r
z = zax[zax < 60]
ne = len(Eax)
nz = len(z)
nc = 3
ns = 2
c = ['HIon', 'LyA', 'Heat']
s = ['y', 'e']

w = np.zeros([ns, nc, nz, ne])
for eid in np.arange(0, ne):
    for zid in np.arange(0, nz):
        for sid in [0, 1]:
            for cid in [0, 1, 2]:
                w_ = Compute_Weigh(
                    idxE=eid, 
                    idxz=Find_idx(z=z[zid]), 
                    spec=s[sid],
                    channel=c[cid])
                w[sid, cid, zid, eid] = w_

x0 = w[Species, 0, :, :]
x1 = w[Species, 1, :, :]
x2 = w[Species, 2, :, :]

plt.rcParams.update({'font.family':'Times'})
plt.rcParams['text.usetex'] = True
fig, axs = plt.subplots(1, 3, sharex = False, sharey = False)
fig.set_size_inches(9, 2.5)
Ek, z = np.meshgrid(Eax, z)

for idx in [0, 1, 2]:
    weigh = w[Species, idx, :, :]
    print(np.shape(weigh))
    c10 = axs[idx].pcolor(Ek, z, weigh, cmap='jet', vmin = 0, vmax = 1, rasterized=True)
    cbar = fig.colorbar(c10, ax=axs[idx])
    axs[idx].tick_params(axis='both', which='both', labelsize = FontSize)
    axs[idx].set_xlabel('Energy [GeV]',fontsize=FontSize,fontname='Times New Roman')
    axs[idx].set_ylabel('Redshift',fontsize=FontSize,fontname='Times New Roman')
    axs[idx].tick_params(axis='both', which='both', labelsize = FontSize)
    #axs[idx].set_xticks(np.linspace(10, 60, 6))
    axs[idx].set_xscale('log')
    #axs[idx].set_yscale('log')
    #axs[idx].set_xlim(11, 50)
    #axs[idx].set_ylim(0.05, 3)
    axs[idx].set_title(c[idx],fontsize=FontSize)
plt.tight_layout()
plt.savefig('/Users/cangtao/Desktop/tmp.pdf')
