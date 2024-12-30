LineWidth = 2
FontSize = 15

# Plot Boost factor
from PyLab import np, plt
d = np.load('../data/BoostTemplates.npz')
z = d['z']
B = d['B']
HMF = d['HMF']
PS = d['PS']
print('PS =',PS)
print('HMF =',HMF)

plt.rcParams.update({
    'text.usetex': True,
    'font.family':'Times',
    'text.latex.preamble': r'\usepackage{newtxtext,newtxmath}'})
fig, axs = plt.subplots(1, 5, sharex = 1, sharey = 1)
fig.set_size_inches(4*5, 4)

color = ['k', 'r', 'b', 'g', 'c']
legends = ['EH', 'BBKS', 'EFSTATHIOU', 'PEEBLES', 'WHITE']
titles = ['PS', 'ST', 'Watson', 'Watson-z']

titles = ['EH', 'BBKS', 'EFSTATHIOU', 'PEEBLES', 'WHITE']
legends = ['PS', 'ST', 'Watson', 'Watson-z']

for pid in PS:
    axs[pid].grid(True, which='both', linewidth = 0.2)
    for hid in HMF:
        axs[pid].semilogy(z, B[hid, pid, :], color[hid], linewidth = LineWidth, label = legends[hid])
    axs[pid].set_xlabel('$z$',fontsize=FontSize,fontname='Times New Roman')
    axs[pid].set_xlim(4, 90)
    axs[pid].tick_params(axis='both', which='both', labelsize = FontSize)
    axs[pid].set_title('Transfer function = ' + titles[pid],fontsize=FontSize)
    axs[pid].set_xticks(np.linspace(10, 90, 9))
    axs[pid].set_yticks(np.logspace(0, 7, 8))
axs[0].set_ylabel('$B$',fontsize=FontSize,fontname='Times New Roman')
axs[0].legend(fontsize=FontSize, loc = 'upper right')
plt.tight_layout()
plt.savefig('/Users/cangtao/Desktop/tmp.pdf')

