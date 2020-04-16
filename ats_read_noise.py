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


def group_image_pairs(file_list):
    image_pair_lists = list(itertools.combinations(file_list, 2))
    return image_pair_lists


def diff_image_stats(img1, img2, sigma_clip=True):
    diff_img = (img1 - img2).flatten()
    if sigma_clip:
        mean, med, stddev = sigma_clipped_stats(diff_img)
        var = stddev**2./2.
    else:
        mean = np.mean(diff_img)
        med = np.median(diff_img)
        var = np.stddev(diff_img)**2./2.
    return mean, med, var

def rebin_image(img, bin_row, bin_col):
    kernel = np.ones((bin_row, bin_col)).astype(np.int)
    c = convolve(img, kernel, mode='valid')
    return c[::bin_row, ::bin_col]

def get_read_noise(pair_list, by_amp=True, roi=[None,None,None,None], bin_row=1, bin_col=1, npairs=250):
    """
    Use provided file list to make a read noise measurement using the diff image method.

    returns mean, variance, read_noise [DN]
    """

    idxs = np.random.choice(np.arange(0,len(pair_list)), npairs, replace=False)
    pair_list = np.array(pair_list)[idxs]
    if by_amp:
        variances = np.zeros((len(pair_list), NUM_AMPS))*np.nan
        means = np.zeros_like(variances)*np.nan
        read_noises = np.zeros_like(variances)*np.nan
        segment_names = np.zeros(NUM_AMPS).astype(np.str)
    else:
        variances = np.zeros(len(pair_list))*np.nan
        means = np.zeros_like(variances)*np.nan
        read_noises = np.zeros_like(variances)*np.nan
        segment_names = [None]

    for i,files in enumerate(pair_list):
        print('Starting set {} of {}...'.format(i,len(pair_list)))
        f1, f2 = files
        d1 = atsi.ATSImage(f1)
        d2 = atsi.ATSImage(f2)
        assert(d1.imagetype=='BIAS'), 'IMAGE {} IS NOT A BIAS'.format(f1)
        assert(d2.imagetype=='BIAS'), 'IMAGE {} IS NOT A BIAS'.format(f2)

        if by_amp:
            for j,a in enumerate(zip(d1.amp_names, d2.amp_names)):
                a1, a2 = a
                assert(a1 == a2), 'SEGMENTS DO NOT MATCH!'
                segment_names[j] = a1
                region1 = d1.amp_images[a1][roi[0]:roi[1],roi[2]:roi[3]]
                region2 = d2.amp_images[a2][roi[0]:roi[1],roi[2]:roi[3]] 

                if bin_row != 1 or bin_col != 1:
                    region1 = rebin_image(region1, bin_row, bin_col)
                    region2 = rebin_image(region2, bin_row, bin_col)
                region1 = region1.astype(np.float)
                region2 = region2.astype(np.float)
                di_mean, di_med, di_var = diff_image_stats(region1, region2)
                variances[i,j] = di_var
                means[i,j] = np.mean((region1+region2).flatten()/2.)
                read_noises[i,j] = np.sqrt(variances[i,j])
                print('Read noise, amp {}: {}'.format(a1, read_noises[i,j]))
        else:
            region1 = d1.image[roi[0]:roi[1],roi[2]:roi[3]]
            region2 = d2.image[roi[0]:roi[1],roi[2]:roi[3]] 
            di_mean, di_med, di_var = diff_image_stats(region1, region2)
            variances[i] = di_var
            means[i] = np.mean((region1+region2).flatten()/2.)
            read_noises[i] = means[i]/variances[i]
            print('Read noise: {}'.format(read_noises[i]))

    return means, variances, read_noises, segment_names


#fnums = np.arange(348, 453, 1)
#fnums = np.arange(248,348,1)
fnums = np.arange(17,191).astype(np.int)

fbase = '/lsstdata/offline/teststand/auxTel/L1Archiver/gen2repo/raw/2020-02-17/2020021700'
fend = '-det000.fits'


flist = [fbase+'%03d'%f+fend for f in fnums]

image_pair_list = np.array(group_image_pairs(flist))

means, variances, read_noises, segment_names = get_read_noise(image_pair_list, by_amp=True, roi=ROI, bin_row=BIN_ROW, bin_col=BIN_COL, npairs=10)

f1 = plt.figure()
plt.ylabel('Number')

for i in range(means.shape[-1]):
    plt.figure(f1.number)
    plt.hist(np.sqrt(variances[:,i]), label='{}'.format(segment_names[i][-2:]), histtype='step', range=[7,10], bins=20)

plt.legend(ncol=2)
plt.xlabel('Read Noise [DN]')

plt.show()


if False:
    for i,name in enumerate(segment_names):
        fname = './read_noise_by_amp/{}.dat'.format(name)
        np.savetxt(fname, np.column_stack((means[:,i], variances[:,i], read_noises[:,i])),
                   header='Mean-Signal-[DN], Variance-[DN^2], Read-noise-[e-/DN]')
    

amp_counts = np.zeros((len(flist), 16)).astype(np.float)
for i,f in enumerate(flist):
    print('Working on image {} of {}'.format(i,len(flist)))
    d = atsi.ATSImage(f)
    for j,amp in enumerate(d.amp_names):
        amp_counts[i,j] = np.median(d.amp_images[amp])

x = np.arange(amp_counts.shape[0])
for j, amp in enumerate(d.amp_names):
    plt.figure()
    plt.plot(x, amp_counts[:,j], '-k.')
    plt.title(amp)
plt.show()

