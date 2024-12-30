LineWidth = 2
FontSize = 15

# Plot Boost factor
from PyLab import np, plt
d = np.load('../data/BoostTemplates.npz')
z = d['z']
B = d['B']
HMF = d['HMF']
PS = d['PS']
print(PS)
print(HMF)

plt.rcParams.update({
    'text.usetex': True,
    'font.family':'Times',
    'text.latex.preamble': r'\usepackage{newtxtext,newtxmath}'})
fig, axs = plt.subplots(1, 4, sharex = 1, sharey = 1)
fig.set_size_inches(4*4, 4)

color = ['k', 'r', 'b', 'g', 'c']
legends = ['EH', 'BBKS', 'EFSTATHIOU', 'PEEBLES', 'WHITE']
titles = ['PS', 'ST', 'Watson', 'Watson-z']

for hid in HMF:
    axs[hid].grid(True, which='both', linewidth = 0.2)
    for pid in PS:
        axs[hid].semilogy(z, B[hid, pid, :], color[pid], linewidth = LineWidth, label = legends[pid])
    axs[hid].set_xlabel('$z$',fontsize=FontSize,fontname='Times New Roman')
    axs[hid].set_xlim(4, 90)
    axs[hid].tick_params(axis='both', which='both', labelsize = FontSize)
    axs[hid].set_title('HMF = ' + titles[hid],fontsize=FontSize)
    axs[hid].set_xticks(np.linspace(10, 90, 9))
    axs[hid].set_yticks(np.logspace(0, 7, 8))
axs[0].set_ylabel('$B$',fontsize=FontSize,fontname='Times New Roman')
axs[0].legend(fontsize=FontSize, loc = 'upper right')
plt.tight_layout()
plt.savefig('/Users/cangtao/Desktop/tmp.pdf')

