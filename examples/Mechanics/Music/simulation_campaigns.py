"""Save campaigns parameters : map oar job ids, dates, frequencies ...
"""
import numpy as np
import os

def create_campaign(freq2job, frequencies, filepath, hosts):
    campaign = {}
    nbfiles = len(freq2job)
    filelist = [None, ] * nbfiles
    newfilelist = [None, ] * nbfiles
    durations = [None, ] * nbfiles
    filepaths = [os.path.join(filepath,
                              'F_' + str(frequencies[i]) + '_id_' + str(freq2job[frequencies[i]][0]))
                 for i in  range(nbfiles)]
    for i in range(nbfiles):
        filelist[i] = 'g_862_' + str(frequencies[i]) + '.h5'
        newfilelist[i] = 'converted_' + filelist[i]
        filelist[i] = os.path.join(filepaths[i], filelist[i])
        newfilelist[i] = os.path.join(filepaths[i], newfilelist[i])
        durations[i] = freq2job[frequencies[i]][1]

    campaign['freqs'] = frequencies
    campaign['files'] = filelist
    campaign['files_converted'] = newfilelist
    campaign['durations'] = durations
    campaign['hosts'] = hosts
    return campaign


#filepath = '/nfs_scratch/perignon/Music'
indices = np.arange(3, 16)
freqs = 2 ** indices * 1960
# Simu on Luke, bass guitar with frets, coeff restit=0.9, 14/11/2017
filepath = './results'
# freq2job = {}
# hosts = [42, 44, 43, 43, 43, 43, 43, 43, 38, 39, 40, 43, 44]

# freq2job[freqs[0]] = (3761741, 5.162)
# freq2job[freqs[1]] = (3761740, 7.941)
# freq2job[freqs[2]] = (3741615, 14.285)
# freq2job[freqs[3]] = (3741614, 27.0571)
# freq2job[freqs[4]] = (3741597, 51.0162399)
# freq2job[freqs[5]] = (3741501, 101.440)
# freq2job[freqs[6]] = (3741500, 200.212)
# freq2job[freqs[7]] = (3741499, 408.23563199999995)
# freq2job[freqs[8]] = (3741498, 1486.282779)
# freq2job[freqs[9]] = (3741496, 2384.006483)
# freq2job[freqs[10]] = (3741495, 5889.027507)
# freq2job[freqs[11]] = (3741492, 6329.552029)
# freq2job[freqs[12]] = (3741487, 12658.475738)
# campaign_14112017 = create_campaign(freq2job, freqs, filepath, hosts)


# # Simu on Luke, bass guitar with frets, coeff restit=0., 16/11/2017
# freq2job = {}
# hosts = [42, 42, 43, 43, 43, 42, 42, 44, 39, 39, 38, 42, 44]

# freq2job[freqs[0]] = (3759710, 4.81778)
# freq2job[freqs[1]] = (3759709,  8.03)
# freq2job[freqs[2]] = (3759708, 15.032)
# freq2job[freqs[3]] = (3759706, 26.468)
# freq2job[freqs[4]] = (3759705, 51.393)
# freq2job[freqs[5]] = (3759697, 105.269)
# freq2job[freqs[6]] = (3759694, 216.863795)
# freq2job[freqs[7]] = (3759691, 404.4368)
# freq2job[freqs[8]] = (3759689, 1203.202637)
# freq2job[freqs[9]] = (3759687, 2379.910844)
# freq2job[freqs[10]] = (3759685, 5941.138048)
# freq2job[freqs[11]] = (3759683,  6573.4847439)
# freq2job[freqs[12]] = (3759682, 12964.30463)
# campaign_16112017 = create_campaign(freq2job, freqs, filepath, hosts)


# # Simu on Luke, fretless bass guitar, coeff restit=0.9, 17/11/2017

# #filepath = os.path.join(filepath, 'Fretless')
# filepath = 'results_fretless'
# freq2job = {}
# hosts = [42, 42, 42, 42, 42, 43, 43, 43, 43, 43, 43, 43, 42]

# freq2job[freqs[0]] = (3759853, 166.051)
# freq2job[freqs[1]] = (3759910, 422.15026300)
# freq2job[freqs[2]] = (3759909, 737.40118)
# freq2job[freqs[3]] = (3759908, 1295.944)
# freq2job[freqs[4]] = (3759907, 2288.943)
# freq2job[freqs[5]] = (3759905, 7076.617)
# freq2job[freqs[6]] = (3759904, 14450.497744)
# freq2job[freqs[7]] = (3759900, 26475.26404)
# freq2job[freqs[8]] = (3759898, 51519.2460)
# freq2job[freqs[9]] = (3759895, 97752.360676)
# freq2job[freqs[10]] = (3759891, 170171.4795)
# freq2job[freqs[11]] = (3759890, 291309.0038)
# freq2job[freqs[12]] = (3759888, 512800.743709)
# campaign_17112017 = create_campaign(freq2job, freqs, filepath, hosts)


# # Simu on Luke, fretless bass guitar, coeff restit=0., 17/11/2017
# freq2job = {}
# hosts = [42, 42, 42, 42, 42, 38, 43, 43, 43, 39, 38, 43, 42]

# freq2job[freqs[0]] = (3760057, 289.98123)
# freq2job[freqs[1]] = (3760056, 486.4436)
# freq2job[freqs[2]] = (3760055, 860.43189)
# freq2job[freqs[3]] = (3760054, 1592.9131)
# freq2job[freqs[4]] = (3760053, 2707.151779)
# freq2job[freqs[5]] = (3760052, 7923.982)
# freq2job[freqs[6]] = (3760051, 14676.945764)
# freq2job[freqs[7]] = (3760050, 27786.139167)
# freq2job[freqs[8]] = (3760049, 49804.66623)
# freq2job[freqs[9]] = (3760048, 83200.4357)
# freq2job[freqs[10]] = (3760047, 211982.44246)
# freq2job[freqs[11]] = (3760046, 288908.6019)
# freq2job[freqs[12]] = (3760045, 504178.2291)
# campaign_17112017_2 = create_campaign(freq2job, freqs, filepath, hosts)




# # Simu on Luke, bass guitar, coeff restit=0.9, 24/11/2017
# freq2job = {}
# hosts = [43, ] * 13

# freq2job[freqs[0]] = (3847227, 5.254)
# freq2job[freqs[1]] = (3847226, 8.896)
# freq2job[freqs[2]] = (3847225, 15.73473)
# freq2job[freqs[3]] = (3847224, 29.5580)
# freq2job[freqs[4]] = (3847223, 57.234)
# freq2job[freqs[5]] = (3847222, 113.259)
# freq2job[freqs[6]] = (3847221, 223.137)
# freq2job[freqs[7]] = (3847220, 444.1789)
# freq2job[freqs[8]] = (3847219, 885.376)
# freq2job[freqs[9]] = (3847218, 1737.819)
# freq2job[freqs[10]] = (3847217, 3439.881) 
# freq2job[freqs[11]] = (3847216, 6526.6391)
# freq2job[freqs[12]] = (3847215, 12785.688522)
# campaign_24112017 = create_campaign(freq2job, freqs, filepath, hosts)

# # Simu on Luke, bass guitar, coeff restit=0., 24/11/2017
# freq2job = {}
# hosts = [42, ] * 13

# freq2job[freqs[0]] = (3847212, 5.58230)
# freq2job[freqs[1]] = (3847211, 9.19436) 
# freq2job[freqs[2]] = (3847210, 16.5804)
# freq2job[freqs[3]] = (3847209, 30.95370999999)
# freq2job[freqs[4]] = (3847208, 58.372619)
# freq2job[freqs[5]] = (3847207, 112.278427)
# freq2job[freqs[6]] = (3847206, 220.387868)
# freq2job[freqs[7]] = (3847205, 473.728463)
# freq2job[freqs[8]] = (3847204, 953.573802)
# freq2job[freqs[9]] = (3847203, 1889.834139)
# freq2job[freqs[10]] = (3847202, 3780.25403)
# freq2job[freqs[11]] = (3847201, 7537.60936)
# freq2job[freqs[12]] = (3847200, 15380.9608)
# campaign_24112017_2 = create_campaign(freq2job, freqs, filepath, hosts)

# # Simu on Luke, fretless bass guitar, coeff restit=0.9 , 27/11/2017
# freq2job = {}
# hosts = [42, ] * 13

# freq2job[freqs[0]] = (3847165, 359.2673)
# freq2job[freqs[1]] = (3847164, 619.7550)
# freq2job[freqs[2]] = (3847163, 1154.076)
# freq2job[freqs[3]] = (3847162, 2145.8427)
# freq2job[freqs[4]] = (3847161, 4215.99275)
# freq2job[freqs[5]] = (3847160, 7832.1030)
# freq2job[freqs[6]] = (3847159, 15052.19231)
# freq2job[freqs[7]] = (3847158, 27637.76544)
# freq2job[freqs[8]] = (3847157, 50594.106563)
# freq2job[freqs[9]] = (3847156, 89999.06627)
# freq2job[freqs[10]] = (3847155, 158867.62378)
# freq2job[freqs[11]] = (3847154, None) # Running)
# freq2job[freqs[12]] = (3847153, None) # Running)
# campaign_27112017 = create_campaign(freq2job, freqs, filepath, hosts)


# # Simu on Luke, fretless bass guitar, coeff restit=0. , 27/11/2017
# freq2job = {}
# hosts = [44, ] * 13

# freq2job[freqs[0]] = (3847252, 256.35580)
# freq2job[freqs[1]] = (3847251, 376.1458)
# freq2job[freqs[2]] = (3847250, 604.6691)
# freq2job[freqs[3]] = (3847249, 1068.560)
# freq2job[freqs[4]] = (3847248, 3783.48298)
# freq2job[freqs[5]] = (3847247, 7267.7895)
# freq2job[freqs[6]] = (3847246, 13774.74347)
# freq2job[freqs[7]] = (3847245, 25960.725063)
# freq2job[freqs[8]] = (3847244, 48870.166304)
# freq2job[freqs[9]] = (3847243, 139390.051575)
# freq2job[freqs[10]] = (3847242, None) # Running)
# freq2job[freqs[11]] = (3847241, None) # Running)
# freq2job[freqs[12]] = (3847240, None) # Running)

# campaign_27112017_2 = create_campaign(freq2job, freqs, filepath, hosts)




# Luke, bass guitar, coeff restit = 0.9, start on 28/11/2017
filepath = './results_bass/2017_11_28'
freq2job = {}
hosts = [42, ] * 13

freq2job[freqs[12]] = (3875402, 123919.2328)
freq2job[freqs[11]] = (3875403, 69255.26833)
freq2job[freqs[10]] = (3875404, 37453.121375)
freq2job[freqs[9]] = (3875405, 19193.3257)
freq2job[freqs[8]] = (3875406, 9366.6517)
freq2job[freqs[7]] = (3875407, 4675.471347)
freq2job[freqs[6]] = (3875408, 2279.230298)
freq2job[freqs[5]] = (3875409, 1079.76208)
freq2job[freqs[4]] = (3875410, 559.481385)
freq2job[freqs[3]] = (3875411, 235.766678)
freq2job[freqs[2]] = (3875412, 83.5007)
freq2job[freqs[1]] = (3875413, 34.7391)
freq2job[freqs[0]] = (3875414, 12.0917)

g2017_11_28_e09 = create_campaign(freq2job, freqs, filepath, hosts)

# Luke, bass guitar, coeff restit = 0., start on 28/11/2017
filepath = './results_bass/2017_11_28'
freq2job = {}
hosts = [44, ] * 13

freq2job[freqs[12]] = (3875417, 120588.15553)
freq2job[freqs[11]] = (3875418, 67350.76821)
freq2job[freqs[10]] = (3875419, 35742.55494)
freq2job[freqs[9]] = (3875420, 18155.32000)
freq2job[freqs[8]] = (3875421, 9054.58146)
freq2job[freqs[7]] = (3875422, 4595.415218)
freq2job[freqs[6]] = (3875423, 2200.718668)
freq2job[freqs[5]] = (3875424, 925.127147)
freq2job[freqs[4]] = (3875425, 456.2572)
freq2job[freqs[3]] = (3875426, 212.18840)
freq2job[freqs[2]] = (3875427, 78.22413)
freq2job[freqs[1]] = (3875428, 27.3738)
freq2job[freqs[0]] = (3875429, 10.4868)

g2017_11_28_e0 = create_campaign(freq2job, freqs, filepath, hosts)

# Luke, bass guitar, coeff restit = 1., start on 28/11/2017
filepath = './results_bass/2017_11_28'
freq2job = {}
hosts = [43, ] * 13

freq2job[freqs[12]] = (3875484, 105420.163)
freq2job[freqs[11]] = (3875485, 60064.23483)
freq2job[freqs[10]] = (3875486, 33331.73031999)
freq2job[freqs[9]] = (3875487, 17431.2319)
freq2job[freqs[8]] = (3875488, 8807.0394)
freq2job[freqs[7]] = (3875489, 4234.7952)
freq2job[freqs[6]] = (3875490, 2008.610674)
freq2job[freqs[5]] = (3875491, 934.37739)
freq2job[freqs[4]] = (3875492, 491.5804)
freq2job[freqs[3]] = (3875493, 214.3197)
freq2job[freqs[2]] = (3875494, 83.6265)
freq2job[freqs[1]] = (3875495, 31.1531)
freq2job[freqs[0]] = (3875496, 5.9595)

g2017_11_28_e1 = create_campaign(freq2job, freqs, filepath, hosts)

                      
# Luke, fretless bass guitar, coeff restit = 0.9, start on 28/11/2017
filepath = './results_fretless/2017_11_28'
freq2job = {}
hosts = [42, ] * 13
freq2job[freqs[12]] = (3876050, 0)
freq2job[freqs[11]] = (3876051, 0)
freq2job[freqs[10]] = (3876052, 0)
freq2job[freqs[9]] = (3876053, 0)
freq2job[freqs[8]] = (3876054, 0)
freq2job[freqs[7]] = (3876055, 0)
freq2job[freqs[6]] = (3876056, 0)
freq2job[freqs[5]] = (3876057, 18871.5758)
freq2job[freqs[4]] = (3876058, 4772.3764)
freq2job[freqs[3]] = (3876059, 2995.850157)
freq2job[freqs[2]] = (3876060, 1827.995)
freq2job[freqs[1]] = (3876061, 910.0303)
freq2job[freqs[0]] = (3876062, 933.8776)


f2017_11_28_e09 = create_campaign(freq2job, freqs, filepath, hosts)

# Luke, fretless bass guitar, coeff restit = 0.9, start on 28/11/2017
filepath = './results_fretless/2017_11_28'
freq2job = {}
hosts = [44, ] * 13
freq2job[freqs[12]] = (3876063, 0)
freq2job[freqs[11]] = (3876064, 0)
freq2job[freqs[10]] = (3876065, 0)
freq2job[freqs[9]] = (3876066, 0)
freq2job[freqs[8]] = (3876067, 0)
freq2job[freqs[7]] = (3876068, 0)
freq2job[freqs[6]] = (3876069, 0)
freq2job[freqs[5]] = (3876070, 21724.836)
freq2job[freqs[4]] = (3876071, 13180.0367)
freq2job[freqs[3]] = (3876072, 3888.4868)
freq2job[freqs[2]] = (3876073, 1467.9752)
freq2job[freqs[1]] = (3876074, 810.9298)
freq2job[freqs[0]] = (3876075, 507.502916)
f2017_11_28_e0 = create_campaign(freq2job, freqs, filepath, hosts)



# Luke, bass guitar, coeff restit = 0.9, start on 30/11/2017
filepath = './results_bass/2017_11_30'
freq2job = {}
hosts = [38, ] * 13
hosts[0] = 42
freq2job[freqs[12]] = (3878785, 14552.025)
freq2job[freqs[11]] = (3877920,  57632.99008)
freq2job[freqs[10]] = (3877921,  5885.640925)
freq2job[freqs[9]] = (3877922, 16464.855114)
freq2job[freqs[8]] = (3877923, 3140.216699)
freq2job[freqs[7]] = (3877924,  2014.87907)
freq2job[freqs[6]] = (3877925,  445.832628)
freq2job[freqs[5]] = (3877926,  256.352)
freq2job[freqs[4]] = (3877927,  320.92811)
freq2job[freqs[3]] = (3877928,  49.077)
freq2job[freqs[2]] = (3877929,  26.742)
freq2job[freqs[1]] = (3877930,  14.95)
freq2job[freqs[0]] = (3877931,  9.495)

g2017_11_30_e09 = create_campaign(freq2job, freqs, filepath, hosts)

# Luke, bass guitar, coeff restit = 0., start on 30/11/2017
filepath = './results_bass/2017_11_30'
freq2job = {}
hosts = [44, ] * 13

freq2job[freqs[12]] = (3877938, 18424.910)
freq2job[freqs[11]] = (3877939, 7833.075)
freq2job[freqs[10]] = (3877940,  4410.211975)
freq2job[freqs[9]] = (3877941, 1756.556)
freq2job[freqs[8]] = (3877942, 888.700)
freq2job[freqs[7]] = (3877943, 1016.859)
freq2job[freqs[6]] = (3877944, 207.6254)
freq2job[freqs[5]] = (3877945, 183.0436)
freq2job[freqs[4]] = (3877946, 55.3991)
freq2job[freqs[3]] = (3877947, 34.29018)
freq2job[freqs[2]] = (3877948, 15.1676)
freq2job[freqs[1]] = (3877949, 8.398)
freq2job[freqs[0]] = (3877950, 5.116)

g2017_11_30_e0 = create_campaign(freq2job, freqs, filepath, hosts)

# Luke, bass guitar, coeff restit = 1., start on 30/11/2017
filepath = './results_bass/2017_11_30'
freq2job = {}
hosts = [44, ] * 13

freq2job[freqs[12]] = (3878147,  37465.54369)
freq2job[freqs[11]] = (3878148,  7500.7054)
freq2job[freqs[10]] = (3878149,  4133.0277)
freq2job[freqs[9]] = (3878150,  5191.068211)
freq2job[freqs[8]] = (3878151,  985.9552)
freq2job[freqs[7]] = (3878152,  472.29017300)
freq2job[freqs[6]] = (3878153,  233.954)
freq2job[freqs[5]] = (3878154,  124.46)
freq2job[freqs[4]] = (3878158,  91.18)
freq2job[freqs[3]] = (3878159,  30.957)
freq2job[freqs[2]] = (3878160,  18.193)
freq2job[freqs[1]] = (3878161,  8.322)
freq2job[freqs[0]] = (3878162,  5.37)

g2017_11_30_e1 = create_campaign(freq2job, freqs, filepath, hosts)


# Luke, fretless bass guitar, coeff restit = 0.9, start on 30/11/2017
filepath = './results_fretless/2017_01_12'
freq2job = {}
hosts = [42, ] * 13

freq2job[freqs[12]] = (3878240, 123919.2328)
freq2job[freqs[11]] = (3878241, 69255.26833)
freq2job[freqs[10]] = (3878242, 37453.121375)
freq2job[freqs[9]] = (3878243, 19193.3257)
freq2job[freqs[8]] = (3878244, 9366.6517)
freq2job[freqs[7]] = (3878245, 4675.471347)
freq2job[freqs[6]] = (3878246, 2279.230298)
freq2job[freqs[5]] = (3878247, 1079.76208)
freq2job[freqs[4]] = (3878248, 559.481385)
freq2job[freqs[3]] = (3878249, 235.766678)
freq2job[freqs[2]] = (3878250, 83.5007)
freq2job[freqs[1]] = (3878251, 34.7391)
freq2job[freqs[0]] = (3878252, 12.0917)

f2017_12_01_e09 = create_campaign(freq2job, freqs, filepath, hosts)


# Luke, fretless bass guitar, coeff restit = 0.9, start on 30/11/2017
filepath = './results_fretless/2017_01_12'
freq2job = {}
hosts = [43, ] * 13

freq2job[freqs[12]] = (3878260, 123919.2328)
freq2job[freqs[11]] = (3878261, 69255.26833)
freq2job[freqs[10]] = (3878262, 37453.121375)
freq2job[freqs[9]] = (3878263, 19193.3257)
freq2job[freqs[8]] = (3878264, 9366.6517)
freq2job[freqs[7]] = (3878265, 4675.471347)
freq2job[freqs[6]] = (3878266, 2279.230298)
freq2job[freqs[5]] = (3878267, 1079.76208)
freq2job[freqs[4]] = (3878268, 559.481385)
freq2job[freqs[3]] = (3878269, 235.766678)
freq2job[freqs[2]] = (3878270, 83.5007)
freq2job[freqs[1]] = (3878271, 34.7391)
freq2job[freqs[0]] = (3878272, 12.0917)

f2017_12_01_e0 = create_campaign(freq2job, freqs, filepath, hosts)
