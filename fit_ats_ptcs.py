import numpy as np
import glob
import matplotlib.pyplot as plt


FIT_MAX = 100000
LOW_LIGHT_LIMIT = 10000
RUNID = 1

run_id_labels = ['2020-02-21', '2020-02-19', '2020-03-13_1', '2020-03-13_2']
limits = np.array([[-0.5,0.75], [-0.5,0.75], [-0.03,0.03], [-0.03,0.03]])

file_list = sorted(glob.glob('../ptc_by_amp/Segment??_binning_1x1_roi_None*'))
rn_file_base = '../read_noise_by_amp/'

fig, ax = plt.subplots(4,4,sharex=True, sharey=True, figsize=(14,14))
ax = ax.flatten()
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off', which='both', length=0)
plt.grid(False)
plt.xlabel('Signal [DN]')
plt.ylabel('S_ADC [e-/DN]', labelpad=12)
plt.title(run_id_labels[RUNID-1])

k_fig = plt.figure()
plt.xlabel('Signal [DN]')
plt.ylabel('K [e-/DN]')

for i,f in enumerate(file_list):
    segment_name = f.split('/')[-1][:9]
    rn_file = rn_file_base + segment_name + '.dat'

    _,_,_,rn = np.loadtxt(rn_file, unpack=True)
    rn = rn.mean()
    means, di_means, var, K, etime, runid = np.loadtxt(f, unpack=True)
    sel = np.where(runid==RUNID)[0]
    means = means[sel]
    var = var[sel]
    etime = etime[sel]
    
    srt = etime.argsort()
    etime = etime[srt]
    light_level = etime[srt]
    light_gain = light_level / light_level.min()
    K = K[srt]
    shot_read_noise = np.sqrt(var[srt])
    signal = means[srt]
    
    # Follow appendix A of Photon Transfer (Janesick)
    shot_noise = np.sqrt(shot_read_noise**2. - rn**2.)
    #Sadc_low = np.mean(signal[signal<LOW_LIGHT_LIMIT] / shot_noise[signal<LOW_LIGHT_LIMIT]**2.)
    #Nadc = np.mean(signal[signal<LOW_LIGHT_LIMIT] / shot_noise[signal<LOW_LIGHT_LIMIT]**2.)
    #s1 = np.mean(signal[signal<LOW_LIGHT_LIMIT])*Sadc_low
    
    Sadc_low = np.mean(signal[0] / shot_noise[0]**2.)
    Nadc = np.mean(signal[0] / shot_noise[0]**2.)
    s1 = np.mean(signal[0])*Sadc_low
    print('e- in first measurement: {}'.format(s1))

    S = s1*light_gain
    Sadc = S / signal
    non_linearity = (Sadc - Sadc_low)/Sadc_low
    
    y = limits[RUNID-1][1]*0.75
    plt.figure(fig.number)
    ax[i].text(x=10000, y=y, s=segment_name, ha='left', va='top', backgroundcolor='white')
    ax[i].plot(means, non_linearity, 'k.')

    plt.figure(k_fig.number)
    plt.scatter(signal, K/K[0], marker='o', label=segment_name[-2:])


plt.figure(fig.number)
for a in ax:
    for i in range(14,18):
        a.axvline(2**i, color='b', ls='--', alpha=0.5)

    a.set_ylim(limits[RUNID-1])
    a.set_xlim(means.min(), means.max())
    a.grid(True)

plt.savefig(run_id_labels[RUNID-1] + '_sadc_nonlin.png')

plt.figure(k_fig.number)
plt.legend(ncol=4)

plt.show()

