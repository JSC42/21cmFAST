UseBoost=1
HMF=0
spec=1
dep_channel=1
idxe = 30
LineWidth = 2
FontSize = 18

from PyLab import plt, np
# A simple script for plotting fc for different PS
EffData = np.load('/Users/cangtao/FileVault/GitHub/21cmFAST/src/py21cmfast/DarkSide/data/EEE_Table.npz')
Ek = EffData['Ek'][idxe]
z = EffData['z']
fc = EffData['fc'][UseBoost, spec, dep_channel, HMF, :, idxe, :]

f0 = fc[0, :]
f1 = fc[1, :]
f2 = fc[2, :]
f3 = fc[3, :]

plt.rcParams.update({
    'text.usetex': True,
    'font.family':'Times',
    'text.latex.preamble': r'\usepackage{newtxtext,newtxmath}'})
fig, ax = plt.subplots()

plt.loglog(1+z, f0, 'k', linewidth=LineWidth)
plt.loglog(1+z, f1, 'r', linewidth=LineWidth)
plt.loglog(1+z, f2, 'b', linewidth=LineWidth)
plt.loglog(1+z, f3, 'g', linewidth=LineWidth)

plt.xlabel('$1+z$',fontsize=FontSize,fontname='Times New Roman')
plt.ylabel('$f_{\mathrm{c}}$',fontsize=FontSize,fontname='Times New Roman')
plt.xticks(size=FontSize)
plt.yticks(size=FontSize)
plt.tight_layout()
plt.savefig('/Users/cangtao/Desktop/tmp.pdf')
