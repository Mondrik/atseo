import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pft
import ats_image as atsi
import ats_helpers as atsh
from scipy.signal import convolve
import itertools
from astropy.stats import sigma_clipped_stats

NUM_AMPS = 16
BIN_ROW = 1
BIN_COL = 1
ROI = [875,1125,375,625] 
ROI = [None, None, None, None]


def get_ptc(pair_list, by_amp=True, roi=[None,None,None,None], bin_row=1, bin_col=1, read_noise=0.):
    """
    Use provided file list to make a photon transfer curve using the diff image method.
    Assumes exposures of same exposure length are grouped, but do *not* need to specify break points.
    E.g., the exposure times can be of the form [e1,e1,e1,e1,e2,e2,e2,e2,...]
    but can not be mixed [e1,e2,e1,e1,e3,e2,...]


    returns mean, variance, K-values (from Janesick)
    """
    exp_times = np.zeros(pair_list.shape[0])
    if by_amp:
        variances = np.zeros((pair_list.shape[0], NUM_AMPS))*np.nan
        means = np.zeros_like(variances)*np.nan
        K_values = np.zeros_like(variances)*np.nan
        segment_names = np.zeros(NUM_AMPS).astype(np.str)
    else:
        variances = np.zeros(pair_list.shape[0])*np.nan
        means = np.zeros_like(variances)*np.nan
        K_values = np.zeros_like(variances)*np.nan
        segment_names = [None]

    remove_rows = []
    for i,files in enumerate(pair_list):
        f1,f2 = files
        print('Working on imgs {}, {}'.format(f1,f2))
        print('image pair: {} of {}'.format(i+1,pair_list.shape[0]))
        d1 = atsi.ATSImage(f1)
        d2 = atsi.ATSImage(f2)
        if d1.exptime != d2.exptime:
            print('Skipping, Exp times didn\'t match...')
            remove_rows.append(i)
            continue
        exp_times[i] = d1.exptime
        print(d1.exptime, d2.exptime)

        if by_amp:
            for j,a in enumerate(zip(d1.amp_names, d2.amp_names)):
                a1, a2 = a
                assert(a1 == a2), 'SEGMENTS DO NOT MATCH!'
                segment_names[j] = a1
                region1 = d1.amp_images[a1][roi[0]:roi[1],roi[2]:roi[3]]
                region2 = d2.amp_images[a2][roi[0]:roi[1],roi[2]:roi[3]] 

                if bin_row != 1 or bin_col != 1:
                    region1 = atsh.rebin_image(region1, bin_row, bin_col)
                    region2 = atsh.rebin_image(region2, bin_row, bin_col)
                di_mean, di_med, di_var = atsh.diff_image_stats(region1, region2)
                variances[i,j] = di_var - read_noise**2.
                means[i,j] = np.mean((region1+region2).flatten()/2.)
                K_values[i,j] = means[i,j]/variances[i,j]
                print('K_value, amp {}: {}'.format(a1[1], K_values[i,j]))
        else:
            region1 = d1.image[roi[0]:roi[1],roi[2]:roi[3]]
            region2 = d2.image[roi[0]:roi[1],roi[2]:roi[3]] 
            di_mean, di_med, di_var = atsh.diff_image_stats(region1, region2)
            variances[i] = vi_var - read_noise**2.
            means[i] = np.mean((region1+region2).flatten()/2.)
            K_values[i] = means[i]/variances[i]
            print('K-value: {}'.format(K_values[i]))

    means = np.delete(means, remove_rows, axis=0)
    variances = np.delete(variances, remove_rows, axis=0)
    K_values = np.delete(K_values, remove_rows, axis=0)
    exp_times = np.delete(exp_times, remove_rows, axis=0)
    return means, variances, K_values, exp_times, segment_names


fend = '-det000.fits'

#PTC Set 1
fbase = '/lsstdata/offline/teststand/auxTel/L1Archiver/gen2repo/raw/2020-02-21/2020022100'
fnums = np.arange(13,61)
flist1 = np.array([fbase+'%03d'%f+fend for f in fnums])


#PTC Set 2
fbase = '/lsstdata/offline/teststand/auxTel/L1Archiver/gen2repo/raw/2020-02-19/2020021900'
fnums = np.arange(41,129)
flist2 = np.array([fbase+'%03d'%f+fend for f in fnums])

#PTC Set 3
fbase = '/lsstdata/offline/teststand/auxTel/L1Archiver/gen2repo/raw/2020-03-13/2020031300'
fnums = np.arange(13,29)
flist3 = np.array([fbase+'%03d'%f+fend for f in fnums])

#PTC Set 4
fbase = '/lsstdata/offline/teststand/auxTel/L1Archiver/gen2repo/raw/2020-03-13/2020031300'
fnums = np.arange(57,97)
flist4 = np.array([fbase+'%03d'%f+fend for f in fnums])

flist = np.concatenate((flist1,flist2,flist3,flist4))



image_pair_list = np.array(atsh.group_image_pairs(flist, by_next=True))

means, variances, K_values, exp_times, segment_names = get_ptc(image_pair_list, by_amp=True, roi=ROI, bin_row=BIN_ROW, bin_col=BIN_COL)

f1 = plt.figure()
plt.title('BINNING: {}x{}'.format(BIN_ROW,BIN_COL))
plt.ylabel('Variance [DN^2]')
plt.yscale('log')
plt.xscale('log')

f2 = plt.figure()
plt.title('BINNING: {}x{}'.format(BIN_ROW,BIN_COL))
plt.ylabel('K [e-/DN]')
plt.ylim(0.9,1.3)

for i in range(means.shape[-1]):
    plt.figure(f1.number)
    plt.plot(means[:,i], variances[:,i], ls='', marker='.', label='{}'.format(segment_names[i]))
    
    plt.figure(f2.number)
    plt.plot(means[:,i], K_values[:,i], ls='', marker='.', label='{}'.format(segment_names[i]))

for f in [f1,f2]:
    plt.figure(f.number)
    plt.legend(ncol=2)
    plt.xlabel('Signal [DN]')

plt.show()

if True:
    for i,name in enumerate(segment_names):
        fname = './ptc_by_amp/{}_binning_{}x{}_roi_{}_{}_{}_{}.dat'.format(
                    name, BIN_ROW, BIN_COL, ROI[0], ROI[1], ROI[2], ROI[3])
        np.savetxt(fname, np.column_stack((means[:,i], variances[:,i], K_values[:,i], exp_times)),
                   header='Mean-Signal-[DN], Variance-[DN^2], K_value-[e-/DN] Exp_time-[s]')
    
    
    
