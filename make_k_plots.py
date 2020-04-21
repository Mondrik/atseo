import numpy as np
import glob
import matplotlib.pyplot as plt

binning = 1
normed = False

file_list = sorted(glob.glob('../ptc_by_amp/Segment??_binning_{:d}x{:d}_roi_None*'.format(binning, binning)))
rn_file_base = '../read_noise_by_amp/'

#k_fig = plt.figure(figsize=(18,9))
fig, ax = plt.subplots(4,4, sharex=True, sharey=True)
ax = ax.flatten()
plt.xlabel('Signal [DN]')
plt.ylabel('K [e-/DN]')

for i,f in enumerate(file_list):
    segment_name = f.split('/')[-1][:9]
    rn_file = rn_file_base + segment_name + '.dat'

    _,_,_,rn = np.loadtxt(rn_file, unpack=True)
    rn = rn.mean()
    print(f'{segment_name} RN: {rn} DN')
    means, di_means, var, K, etime, runid = np.loadtxt(f, unpack=True)
    
    srt = etime.argsort()
    etime = etime[srt]
    K = K[srt]
    signal = means[srt]
    K2 = signal / (var[srt] - (binning*rn)**2.)
    #plt.figure(k_fig.number)
    ax[i].plot(etime, signal, marker='o', ls='')
#    if normed:
#        plt.scatter(signal, K2/K2[0], marker='o', label=segment_name[-2:])
#    else:
#        plt.scatter(signal, K2, marker='o', label=segment_name[-2:])
#
#
#for i in [14,15,16]:
#    plt.axvline(binning*binning*2**i-1, ls='--', color='k')
#
#plt.ylim(1, 1.35)
#plt.title('{}x{} Binning'.format(binning, binning))
#plt.figure(k_fig.number)
#plt.legend(ncol=4)
#plt.tight_layout()
#
#if normed:
#    fname = 'K_values_{}x{}_normed.png'.format(binning, binning)
#else:
#    fname = 'K_values_{}x{}.png'.format(binning, binning)
#
#plt.savefig(fname)
plt.show()
#
