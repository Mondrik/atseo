import numpy as np
import glob
import matplotlib.pyplot as plt


FIT_MAX = 110000
LOW_LIGHT_LIMIT = 10000

file_list = glob.glob('./ptc_by_amp/Segment??_binning_1x1_roi_None*')
rn_file_base = './read_noise_by_amp/'

for f in file_list:
    segment_name = f.split('/')[-1][:9]
    rn_file = rn_file_base + segment_name + '.dat'

    _,_,rn = np.loadtxt(rn_file, unpack=True)
    rn = rn.mean()
    mean, var, K, etime = np.loadtxt(f, unpack=True)
    
    srt = etime.argsort()
    light_level = etime[srt]
    light_gain = light_level / light_level.min()
    K = K[srt]
    shot_read_noise = np.sqrt(var[srt])
    signal = mean[srt]
    
    # Follow appendix A of Photon Transfer (Janesick)
    shot_noise = np.sqrt(shot_read_noise**2. - rn**2.)
    Sadc_low = np.mean(signal[signal<LOW_LIGHT_LIMIT] / shot_noise[signal<LOW_LIGHT_LIMIT]**2.)
    Nadc = np.mean(signal[signal<LOW_LIGHT_LIMIT] / shot_noise[signal<LOW_LIGHT_LIMIT]**2.)
    s1 = np.mean(signal[signal<LOW_LIGHT_LIMIT])*Sadc_low
    print('e- in first measurement: {}'.format(s1))

    S = s1*light_gain
    Sadc = S / signal
    non_linearity = (Sadc - Sadc_low)/Sadc_low
    
    plt.figure()
    plt.title(segment_name)
    plt.plot(mean, non_linearity, 'k.')
    plt.show(block=False)
