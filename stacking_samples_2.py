# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 09:46:02 2015

@author: AmberHornsby
"""

import numpy as np
import glob
import matplotlib.pyplot as plt
#functions#####################################################################
def stack(empty1,empty2,dir1,data,name,fraction,ra_col,dec_col):
    #leave these blank so initialises empty array 
    empty1 = []
    empty2 = []
    # fill in for ra and dec however you splice your array
    for n in range(len(data)):
        empty1.append('samples_*_'+str(data[n,ra_col])+'_'+str(data[n,dec_col])+'_*.npy')
        empty2.append('lnprob_run_*_'+str(data[n,ra_col])+'_'+str(data[n,dec_col])+'_*.npy')
    
    print 'number of galaxies', 
    print len(data)
    
    #define your bins and initialise empty stacked arrays
    X = np.linspace(0, 14, 100)
    Y = np.linspace(0, 4, 100)
    sum1 =np.zeros((len(X)-1, len(Y)-1))
    count = 0
    
    bestfit = np.zeros((1,6))
    
    for j in range(len(empty1)):
        print name, (j/float(len(empty1)))*100, '% complete'
        try:
            s = np.load(glob.glob(dir1+empty1[j])[0])
            p = np.exp(np.load(glob.glob(dir1+empty2[j])[0]))
            s = s[np.where(p>0.2)]
            p = p[np.where(p>0.2)]
            if len(s) == 0:
               pass
            else:
                bestfit = np.append(bestfit, np.array(map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(s, [16,50,84],axis=0)))).reshape(1,6), axis=0)
                Hs, Xs, Ys = np.histogram2d(s[:,0], s[:,1], bins=(X, Y), normed=True, weights=np.log(p))
                #change disc and smooth to whichever wanted to weight by
                sum1 += Hs*fraction[j]
    
        except IndexError:
            count+=1
            print count
            pass
    
    #fill in save names
    np.save('best_fit_'+ name+'.npy', bestfit)
    np.save('stack_' +name+'.npy', sum1)
    return sum1
#directories where bulgless and bulged stored##################################
dir1='/Users/AmberHornsby/github/starpy/Bulgeless_galaxies/'
dir2 = '/Users/AmberHornsby/github/starpy/bulged_gal/'
#Table where data is stored####################################################
data_less = np.genfromtxt('bulgeless_need.txt', delimiter = ',')
data_bulged = np.genfromtxt('bulged_need.txt', delimiter = ',')
fraction_less = np.genfromtxt('bulgeless_fraction.txt', delimiter = ',')
fraction_bulged = np.genfromtxt('bulged_fraction.txt', delimiter = ',')
#stacking samples##############################################################
#st_less = []
#ln_p_less = []
#sum_less = stack(st_less,ln_p_less,dir1,data_less,str('bulgeless'),fraction_less,6,7)
#st_bulged = []
#ln_p_bulged = []
#sum_bulged = stack(st_bulged,ln_p_bulged,dir2,data_bulged,str('bulged'),fraction_bulged,6,7)
plt.figure()
plt.subplot(221)
plt.imshow(sum_bulged.T, origin = 'lower')
plt.colorbar()
plt.ylabel("Bulged")
plt.subplot(222)
plt.imshow(sum_less.T, origin = 'lower')
plt.colorbar()
plt.ylabel("Bulgeless")
plt.subplot(223)
plt.ylabel("log(Bulged)")
plt.imshow(np.log(sum_bulged.T), origin = 'lower')
plt.colorbar()
plt.subplot(224)
plt.ylabel("log(Bulgeless)")
plt.imshow(np.log(sum_less.T), origin = 'lower')
plt.colorbar()

