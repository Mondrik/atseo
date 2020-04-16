import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pft
import ats_image as atsh
from scipy.signal import convolve
import itertools


NUM_AMPS = 16
BIN_ROW = 1
BIN_COL = 1
ROI = [None, None, None, None] 


def group_image_pairs(file_list):
    image_pair_lists = list(itertools.combinations(file_list, 2))
    return image_pair_lists


def diff_image_variance(img1, img2):
    diff_img = img1 - img2
    var = np.std(diff_img.flatten())**2./2.
    return var


def rebin_image(img, bin_row, bin_col):
    kernel = np.ones((bin_row, bin_col)).astype(np.int)
    c = convolve(img, kernel, mode='valid')
    return c[::bin_row, ::bin_col]


#fnums = np.arange(348, 453, 1)
#fnums = np.arange(248,348,1)
fnums = np.arange(17,191).astype(np.int)

fbase = '/lsstdata/offline/teststand/auxTel/L1Archiver/gen2repo/raw/2020-02-17/2020021700'
fend = '-det000.fits'


flist = [fbase+'%03d'%f+fend for f in fnums]

image_pair_list = np.array(group_image_pairs(flist))

overscan_means = np.zeros((len(flist), 16))
bias_means = np.zeros_like(overscan_means)
obs_times = np.zeros(len(flist))

cm = plt.cm.get_cmap('jet', 16)
fig = plt.figure()

for i,f in enumerate(flist):
    print('Working on {} of {}'.format(i,len(flist)))
    d = atsh.ATSImage(f)
    obs_times[i] = d.dateobs.mjd
    for a,amp in enumerate(d.amp_names):
        overscan_means[i,a] = d.amp_overscan_means[amp]
        bias_means[i,a] = d.amp_images[amp].flatten().mean()
        if amp == 'Segment11':
            plt.imshow(d.image, origin='lower', vmin=-50, vmax=50)
            plt.show(block=True)


for i in range(overscan_means.shape[-1]):
        plt.plot(obs_times - obs_times.min(), overscan_means[:,i]-overscan_means[:,i].mean(), marker='.', color=cm(i), label=d.amp_names[i])

plt.legend(ncol=4)
plt.xlabel('Obs Time [d]')
plt.ylabel('Overscan Residual')

fig2 = plt.figure()
for i in range(bias_means.shape[-1]):
        plt.plot(obs_times - obs_times.min(), bias_means[:,i]-bias_means[:,i].mean(), marker='.', color=cm(i), label=d.amp_names[i])

plt.legend(ncol=4)
plt.xlabel('Obs Time [d]')
plt.ylabel('Bias Residual')
plt.show()


plt.show()

