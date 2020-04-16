import astropy.io.fits as pft
import numpy as np
from astropy.time import Time
import itertools
from astropy.stats import sigma_clipped_stats
from scipy.signal import convolve


def group_image_pairs(file_list, by_next=False, by_exptime=False):
    if by_next:
        image_pair_lists = [[file_list[2*i], file_list[2*i+1]] for i in range(int(len(file_list)/2)-1)]
    elif by_exptime:
        exp_times = np.array([np.float(pft.open(f)[0].header['EXPTIME']) for f in flist])
        print(list(set(exp_times)))
        image_pair_lists = []
        for e in exp_times:
            i = np.where(exp_times == e)[0]
            imlist = file_list[i]
            combinations = list(itertools.combinations(imlist,2))
            image_pair_lists += combinations
    else:
        image_pair_lists = list(itertools.combinations(file_list,2))
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


