#%%
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import numpy as np
from multicarrier_fitting import HRMR_PPMS
# %%
plt.rcParams['figure.subplot.bottom'] = 0.1
plt.rcParams['figure.subplot.left'] = 0.1
plt.rcParams["font.size"] = 8.5
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = 'stix'
# %%
files = [
    "dummy_data/MR_HR_5K_500uA.dat",
    "dummy_data/MR_HR_100K_500uA.dat"
]
labels = ["5K", "100K" ]
# %%
popt_list = []
pcov_list = []

for file, label in zip(files, labels):
    data = HRMR_PPMS(file)
    data.remove_even()
    data.to_resistivity(thickness_nm=28.5)
    h_, sigma_xx, sigma_yx = data._sigma()
    popt, pcov = HRMR_PPMS.multi_carrier_fit(h_, sigma_xx, sigma_yx , p0=[1e26, 1e-2, -1e26, -1e-2])
    popt_list.append(popt)
    pcov_list.append(pcov)
    sigma_xx_sigma_yx = HRMR_PPMS.multi_carrier_simu(np.concatenate((h_, h_)), *popt)/1e4

    fig = plt.figure(figsize=(2.3,3.2), dpi=300)
    plus = lambda x, y: x + y
    ax1 = fig.add_subplot(2,1,1)
    ax1.scatter(h_, (sigma_xx)/1e4, label="data", marker=".", s=30, color="blue")
    ax1.plot(h_, sigma_xx_sigma_yx[:len(sigma_xx_sigma_yx)//2], label="2-carrier fit", color="black")
    ax1.tick_params(labelbottom=False)
    ax1.set_ylabel(r"$\sigma_{xx}\ \mathrm{(10^4 S/m)}$")
    ax1.legend(handlelength=0.8, handletextpad=0.2)

    ax2 = fig.add_subplot(2,1,2)
    ax2.scatter(h_, (sigma_yx)/1e4, label="data", marker=".", s=30, color="blue")
    ax2.plot(h_, sigma_xx_sigma_yx[len(sigma_xx_sigma_yx)//2:], label="2-carrier fit", color='black')
    ax2.set_ylabel(r"$\sigma_{yx}\ \mathrm{(10^4 S/m)}$")
    ax2.set_xlabel(r"$μ_0H$ (T)")
    ax2.legend(handlelength=0.8, handletextpad=0.2)
    plt.tight_layout(pad=0.2)
    plt.show()
# %%
popt_array = np.stack(popt_list).T
temp_list = [5, 100]
err = np.array([np.sqrt(np.diag(pcov)) for pcov in pcov_list]).T

fig = plt.figure(figsize=(2.3,3.2), dpi=300)
def plot_all(ax: plt.Axes, x, y, cov, color, label):
    ax.errorbar(x, y, yerr = cov, capsize=2, fmt='.', markersize=5, ecolor=color, markeredgecolor = color, color=color)
    ax.plot(x, y, label=label, color=color)

ax1 = fig.add_subplot(2,1,1)
plot_all(ax1, temp_list, popt_array[0]*1e-26, err[0]*1e-26, 'tab:blue', "n-type")
plot_all(ax1, temp_list, -popt_array[2]*1e-26, err[2]*1e-26, 'tab:red', "p-type")
ax1.legend(frameon=False, bbox_to_anchor=(1, 1), loc='upper right', handlelength=0.8, handletextpad=0.2)
ax1.set_ylabel(r"Carrier density $\mathrm{(10^{20} /cm^3)}$")
plt.setp(ax1.get_xticklabels(), visible=False)

ax2 = fig.add_subplot(2,1,2)
plot_all(ax2, temp_list, popt_array[1]*1e4, err[1]*1e4, 'tab:blue', "n-type")
plot_all(ax2, temp_list, -popt_array[3]*1e4, err[3]*1e4, 'tab:red', "p-type")
ax2.set_xlabel("Temperature (K)")
ax2.set_ylabel(r"Mobility $\mathrm{(cm^2/V^{-1}s^{-1})}$")
ax2.legend(frameon=False, bbox_to_anchor=(1, 1), loc='upper right', handlelength=0.8, handletextpad=0.2)

plt.tight_layout(pad=0.2)
plt.subplots_adjust(hspace=.0)
plt.show()
# plt.savefig("carrier.eps")
# %%
usercmap = plt.get_cmap('gist_rainbow')
cNorm  = colors.Normalize(vmin=0, vmax=6)
scalarMap = cm.ScalarMappable(norm=cNorm, cmap=usercmap)
# %% MRのプロット
plt.figure(figsize=(1.8,2.1), dpi=300)
j=0
for file, label in zip(files, labels):
    data = HRMR_PPMS(file)
    data.to_resistivity(28.5)
    mr = data.mr()
    plt.plot(data.data[0], mr, label=label, color=scalarMap.to_rgba(j), marker=".", markersize=2)
    j += 1
plt.ylabel("MR (%)")
plt.xlabel(r"$μ_0H$ (T)")
plt.legend(frameon=False, handlelength=0.8, handletextpad=0.2, ncol=2, columnspacing=0.2)
plt.tight_layout(pad=0.2)
plt.subplots_adjust(left=0.26)
# plt.savefig("MR.eps")
plt.show()
# %% HRのプロット
plt.figure(figsize=(1.8,2.1), dpi=300)
j=0
for file, label in zip(files, labels):
    data = HRMR_PPMS(file)
    data.remove_even()
    data.to_resistivity(28.5)
    plt.plot(data.data[0], data.data[2]*1e2, label=label, color=scalarMap.to_rgba(j), marker=".", markersize=2)
    j += 1
ax1.tick_params(labelbottom=False)
ax1.set_ylabel(r"$\rho_{xx}\ \mathrm{(10^4 S/m)}$")
ax1.legend(handlelength=0.8, handletextpad=0.2)
ax2.set_ylabel(r"$\rho_{yx}\ \mathrm{(10^4 S/m)}$")
ax2.set_xlabel(r"$μ_0H$ (T)")
plt.ylabel(r"$ρ_{yx}$ (μΩ cm)")
plt.xlabel(r"$μ_0H$ (T)")
#plt.yticks([-0.2,-0.1,0,0.1,0.2],[-0.2,-0.1,0,0.1,0.2])
plt.tight_layout(pad=0.2)
plt.legend(frameon=False, handlelength=0.8, handletextpad=0.2, ncol=2, columnspacing=0.2)
# plt.savefig("HR.eps")
plt.show()
# %%
