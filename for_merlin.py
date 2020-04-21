import numpy as np

fend = '-det000.fits'

#PTC Set 1
fbase = '/lsstdata/offline/teststand/auxTel/L1Archiver/gen2repo/raw/2020-02-21/2020022100'
fnums = np.arange(13,61)
ptc_run_1 = np.ones_like(fnums).astype(int)*1.
flist1 = np.array([fbase+'%03d'%f+fend for f in fnums])


#PTC Set 2
fbase = '/lsstdata/offline/teststand/auxTel/L1Archiver/gen2repo/raw/2020-02-19/2020021900'
fnums = np.arange(41,129)
ptc_run_2 = np.ones_like(fnums).astype(int)*2
flist2 = np.array([fbase+'%03d'%f+fend for f in fnums])

#PTC Set 3
fbase = '/lsstdata/offline/teststand/auxTel/L1Archiver/gen2repo/raw/2020-03-13/2020031300'
fnums = np.arange(13,29)
ptc_run_3 = np.ones_like(fnums).astype(int)*3
flist3 = np.array([fbase+'%03d'%f+fend for f in fnums])

#PTC Set 4
fbase = '/lsstdata/offline/teststand/auxTel/L1Archiver/gen2repo/raw/2020-03-13/2020031300'
fnums = np.arange(57,97)
ptc_run_4 = np.ones_like(fnums).astype(int)*4
flist4 = np.array([fbase+'%03d'%f+fend for f in fnums])

flist = np.concatenate((flist1,flist2,flist3,flist4))

print(flist)
