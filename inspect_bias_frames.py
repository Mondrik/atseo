import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pft
import ats_image as atsi
import ats_helpers as atsh
from scipy.signal import convolve
import itertools


NUM_AMPS = 16
BIN_ROW = 1
BIN_COL = 1
ROI = [None, None, None, None] 


#fnums = np.arange(348, 453, 1)
#fnums = np.arange(248,348,1)
fnums = np.arange(17,191).astype(np.int)

fbase = '/lsstdata/offline/teststand/auxTel/L1Archiver/gen2repo/raw/2020-02-17/2020021700'
fend = '-det000.fits'


flist = [fbase+'%03d'%f+fend for f in fnums]

image_pair_list = np.array(atsh.group_image_pairs(flist))

overscan_means = np.zeros((len(flist), 16))
bias_means = np.zeros_like(overscan_means)
obs_times = np.zeros(len(flist))

cm = plt.cm.get_cmap('jet', 16)
fig = plt.figure()

for i,f in enumerate(flist):
    print('Working on {} of {}'.format(i,len(flist)))
    d = atsi.ATSImage(f)
    obs_times[i] = d.dateobs.mjd
    for a,amp in enumerate(d.amp_names):
        overscan_means[i,a] = d.amp_overscan_means[amp]
        bias_means[i,a] = d.amp_images[amp].flatten().mean()
#        if amp == 'Segment11' and False:
#            plt.imshow(d.image, origin='lower', vmin=-50, vmax=50)
#            plt.show(block=True)
#        plt.hist(d.amp_images[amp].flatten(), range=[-50,50], bins=100)
#    plt.show(block=True)


for i in range(overscan_means.shape[-1]):
    plt.plot(obs_times - obs_times.min(), overscan_means[:,i]-overscan_means[:,i].mean(), marker='.', color=cm(i), label=d.amp_names[i])
    plt.plot(obs_times - obs_times.min(), overscan_means[:,i]-overscan_means[:,i].mean(), marker='.', color=cm(i), label=d.amp_names[i])


plt.legend(ncol=4)
plt.xlabel('Obs Time [d]')
plt.ylabel('Overscan Residual')

fig2 = plt.figure()
for i in range(bias_means.shape[-1]):
    #plt.plot(obs_times - obs_times.min(), bias_means[:,i]-bias_means[:,i].mean(), marker='.', color=cm(i), label=d.amp_names[i])
    plt.plot(bias_means[:,i]-bias_means[:,i].mean(), marker='.', color=cm(i), label=d.amp_names[i])


plt.legend(ncol=4)
plt.xlabel('Obs Time [d]')
plt.ylabel('Bias Residual')
plt.show()


plt.show()

