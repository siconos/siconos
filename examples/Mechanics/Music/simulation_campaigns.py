"""Save campaigns parameters : map oar job ids, dates, frequencies ...
"""
import numpy as np
import os

filepath = '/nfs_scratch/perignon/Music'
indices = np.arange(3, 16)
freqs = 2 ** indices * 1960
# Simu on Luke, bass guitar with frets, coeff restit=0.9, 14/11/2017
campaign_14112017 = {}
freq2job = {}

freq2job[freqs[0]] = 3761741
freq2job[freqs[1]] = 3761740
freq2job[freqs[2]] = 3741615
freq2job[freqs[3]] = 3741614
freq2job[freqs[4]] = 3741597
freq2job[freqs[5]] = 3741501
freq2job[freqs[6]] = 3741500
freq2job[freqs[7]] = 3741499
freq2job[freqs[8]] = 3741498
freq2job[freqs[9]] = 3741496
freq2job[freqs[10]] = 3741495
freq2job[freqs[11]] = 3741492
freq2job[freqs[12]] = 3741487
nbfiles = len(freq2job)
filelist = ['F_' + str(freqs[i]) + '_id_' + str(freq2job[freqs[i]]) for i in  range(nbfiles)]
for i in range(nbfiles):
    filelist[i] = os.path.join(filelist[i],'g_862_' + str(freqs[i]) + '.h5')
    filelist[i] = os.path.join(filepath, filelist[i])
campaign_14112017['freqs'] = freqs
campaign_14112017['files'] = filelist



# Simu on Luke, bass guitar with frets, coeff restit=0., 16/11/2017
campaign_16112017 = {}
freq2job = {}

freq2job[freqs[0]] = 3759710
freq2job[freqs[1]] = 3759709
freq2job[freqs[2]] = 3759708
freq2job[freqs[3]] = 3759706
freq2job[freqs[4]] = 3759705
freq2job[freqs[5]] = 3759697
freq2job[freqs[6]] = 3759694
freq2job[freqs[7]] = 3759691
freq2job[freqs[8]] = 3759689
freq2job[freqs[9]] = 3759687
freq2job[freqs[10]] = 3759685
freq2job[freqs[11]] = 3759683
freq2job[freqs[12]] = 3759682
filelist = ['F_' + str(freqs[i]) + '_id_' + str(freq2job[freqs[i]]) for i in  range(nbfiles)]
for i in range(nbfiles):
    filelist[i] = os.path.join(filelist[i],'g_862_' + str(freqs[i]) + '.h5')
    filelist[i] = os.path.join(filepath, filelist[i])
campaign_16112017['freqs'] = freqs
campaign_16112017['files'] = filelist



# Simu on Luke, fretless bass guitar, coeff restit=0.9, 17/11/2017
filepath = os.path.join(filepath, 'Fretless')

campaign_17112017 = {}
freq2job = {}

freq2job[freqs[0]] = 3759853
freq2job[freqs[1]] = 3759910
freq2job[freqs[2]] = 3759909
freq2job[freqs[3]] = 3759908
freq2job[freqs[4]] = 3759907
freq2job[freqs[5]] = 3759905
freq2job[freqs[6]] = 3759904
freq2job[freqs[7]] = 3759900
freq2job[freqs[8]] = 3759898
freq2job[freqs[9]] = 3759895
freq2job[freqs[10]] = 3759891
freq2job[freqs[11]] = 3759890
freq2job[freqs[12]] = None #3759888
filelist = ['F_' + str(freqs[i]) + '_id_' + str(freq2job[freqs[i]]) for i in  range(nbfiles)]
for i in range(nbfiles):
    filelist[i] = os.path.join(filelist[i],'g_862_' + str(freqs[i]) + '.h5')
    filelist[i] = os.path.join(filepath, filelist[i])
campaign_17112017['freqs'] = freqs
campaign_17112017['files'] = filelist


# Simu on Luke, fretless bass guitar, coeff restit=0., 17/11/2017
campaign_17112017_2 = {}
freq2job = {}

freq2job[freqs[0]] = 3760057
freq2job[freqs[1]] = 3760056
freq2job[freqs[2]] = 3760055
freq2job[freqs[3]] = 3760054
freq2job[freqs[4]] = 3760053
freq2job[freqs[5]] = 3760052
freq2job[freqs[6]] = 3760051
freq2job[freqs[7]] = 3760050
freq2job[freqs[8]] = 3760049
freq2job[freqs[9]] = 3760048
freq2job[freqs[10]] = 3760047
freq2job[freqs[11]] = 3760046
freq2job[freqs[12]] = None # 3760045
filelist = ['F_' + str(freqs[i]) + '_id_' + str(freq2job[freqs[i]]) for i in  range(nbfiles)]
for i in range(nbfiles):
    filelist[i] = os.path.join(filelist[i],'g_862_' + str(freqs[i]) + '.h5')
    filelist[i] = os.path.join(filepath, filelist[i])
campaign_17112017_2['freqs'] = freqs
campaign_17112017_2['files'] = filelist


