# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 14:33:35 2015

@author: AmberHornsby
"""
import matplotlib.pyplot as plt
import numpy as np
import numpy
print numpy.__version__
###############################################################################
#10:10 Classifying galaxies and investigating evolution
###############################################################################
plt.close('all')
#functions#####################################################################
#defining green valley
def green_valley(Mr):
    return 2.06 - 0.244*np.tanh((Mr + 20.07)/1.09)
#defining red galaxies
def red(Mr,colour):
    return colour > green_valley(Mr)
#true false array
def TF_array(R,galaxy,bins_array):
    c = []
    for i in xrange(bins_array.shape[0]):
        #digitize calculates hist bins, in1d identifies where two arrays have 
        #same value
        z = np.where(np.in1d(bins_array[np.digitize(R[galaxy],bins_array)],
                               bins_array[i]))
        c.append(z[0])
    return np.asarray(c)
#Data Input####################################################################
read_data = False
if read_data:
    data = np.genfromtxt('combined.csv', delimiter = ',',skip_header = 1)
    np.save('data', data)
data = np.load('data.npy')
'''
Column 1 = dr8  ;   Column 3 = ra_1     ;   Column 4 = dec
Column 21 = disk    ;   Column 43 = face on     ;   Column 69 = no bulge
Column 75 = just noticeable     ;   Column 81 = obvious     
Column 85 = dominant    ;   Column 278 = redshift   ;   Column 280 = U
Column 281 = G  ;   Column 282 = R  ;   Column 285 = error in U
Column 287 = error in R
'''
np.savetxt("cut.csv", data[0:500], delimiter = ",")
#selecting data ###############################################################
disk = data[...,21]
face_on = data[...,43]
no_bulge = data[...,69]
just_noticeable = data[...,75]
obvious = data[...,81]
dominant = data[...,85]
redshift = data[...,278]
U = data[...,280]
G = data[...,281]
R = data[...,282]
errU = data[...,285]
errR = data[...,287]

#maths#########################################################################
gminusr = G-R
uminusr = U-R
size = 2.5
g = np.linspace(-25,-10,num=1000)
gv = green_valley(g)
#face on cut###################################################################
lim_face_on = 0.430
face_on_loc = np.where(face_on>lim_face_on)
x = R[face_on_loc].shape[0]
y = R.shape[0]
print "total number of galaxies before cuts = %.2d" % y
print "number of face on galaxies = %.2d" % x
#bulge/bulgeless cut###########################################################
lim_bulge = 0.715
#true\false array
no_bulge_cut= np.logical_and(face_on>lim_face_on, no_bulge>lim_bulge)
just_not_cut = np.logical_and(face_on>lim_face_on, just_noticeable>lim_bulge)
obvious_cut = np.logical_and(face_on>lim_face_on, obvious>lim_bulge)
dominant_cut = np.logical_and(face_on>lim_face_on, dominant>lim_bulge)
#test plot bulge/bulgeless galaxies############################################
plt.figure() #figure 1
ax1 = plt.gca()
ax1.invert_xaxis()
plt.plot(R[np.where(no_bulge_cut)], uminusr[np.where(no_bulge_cut)], 'mo', 
           markersize = size,label = "no bulge")
plt.plot(R[np.where(just_not_cut)], uminusr[np.where(just_not_cut)], 'co', 
           markersize = size, label = "just noticeable")
plt.plot(R[np.where(obvious_cut)], uminusr[np.where(obvious_cut)], 'bo', 
           markersize = size, label = "obvious")    
plt.plot(R[np.where(dominant_cut)], uminusr[np.where(dominant_cut)], 'ro', 
           markersize = size, label = "dominant")
plt.plot(g,gv,'g-', linewidth = 3, label = "green valley divider")
plt.xlabel("Absolute R")
plt.ylabel("U - R")
plt.legend(loc="best")
plt.xlim(-16,-25)
plt.ylim(0.5,4)
#number of galaxies############################################################
x = R[np.where(no_bulge_cut)].shape[0]
y = R[np.where(just_not_cut)].shape[0]
z = R[np.where(obvious_cut)].shape[0]
a = R[np.where(dominant_cut)].shape[0]
tot = x + a + z +y
print "number of bulgeless galaxies = %.2d" % x
print "number of just noticeable galaxies = %.2d" % y
print "number of obvious galaxies = %.2d" % z
print "number of dominant galaxies = %.2d" % a
print "total number of galaxies in bulge/bulgeless catagory = %.2d" % tot
#red galaxy cut################################################################
red= red(R,uminusr)
red_bulgeless =  np.where(np.logical_and(red,no_bulge_cut))
red_just_not =  np.where(np.logical_and(red,just_not_cut))
red_obvious =  np.where(np.logical_and(red,obvious_cut))
red_dominant = np.where(np.logical_and(red,dominant_cut))
#adding up all bulged galaxies
red_bulged = np.where(np.logical_or(np.logical_and(red,obvious_cut),
                   np.logical_and(red,dominant_cut)) )    
#saving red galaxy cuts########################################################
redbulgelessrows = data[np.logical_and(red,no_bulge_cut)]   
redbulgedrows = data[np.logical_or(np.logical_and(red,obvious_cut),
                   np.logical_and(red,dominant_cut))]
np.savetxt("red_bulgeless.csv", redbulgelessrows, delimiter = ",")
np.savetxt("red_bulged.csv", redbulgedrows, delimiter = ",")
#test plot red galaxies########################################################
plt.figure() #figure 2
ax = plt.gca()
ax.invert_xaxis()
plt.plot(R[red_bulgeless], uminusr[red_bulgeless], 'mo', markersize = size,
         label = "no bulge")
plt.plot(R[red_just_not], uminusr[red_just_not], 'co', markersize = size,
         label = "just noticeable")
plt.plot(R[red_obvious], uminusr[red_obvious], 'bo', markersize = size,
         label = "obvious")    
plt.plot(R[red_dominant], uminusr[red_dominant], 'ro', markersize = size,
         label = "dominant")
plt.plot(g,gv,'g-', linewidth = 3, label = "green valley divider")
plt.xlabel("Absolute R")
plt.ylabel("U - R")
plt.legend(loc="best")
plt.xlim(-16,-25)
plt.ylim(0.5,4)
#number of red galaxies########################################################
x = R[red_bulgeless].shape[0]
y = R[red_just_not].shape[0]
z = R[red_obvious].shape[0]
a = R[red_dominant].shape[0]
tot = x + a + z + y
print "number of red bulgeless galaxies = %.2d" % x
print "number of red just noticeable galaxies = %.2d" % y
print "number of red obvious galaxies = %.2d" % z
print "number of red dominant galaxies = %.2d" % a
print "total number of galaxies in red bulge/bulgeless catagory = %.2d" % tot
#plotting magnitude histograms#################################################
#importing NUV data from bulged/bulgeless galaxies
#red_bulged = np.genfromtxt("bulged.csv", delimiter = ",", skiprows = 1)
#red_bulgeless = np.genfromtxt("bulgeless.csv", delimiter = ",", skiprows = 1)
###############################################################################
bins = 24
bins_array = np.linspace(-24,-12,num = 13)
plt.figure() #figure 3
plt.subplot(121)
plt.hist(R[red_bulgeless], bins=bins, normed = True, label = "bulgeless")
plt.legend(loc="best")
plt.subplot(122)
plt.hist(R[red_bulged], bins=bins, normed = True, label = "bulged")
plt.legend(loc="best")
#extracting bulgeless data from bins ##########################################
TF_array_bulgeless = TF_array(R,red_bulgeless,bins_array)  
#extracting bulged data from bins #############################################
TF_array_bulged = TF_array(R,red_bulged, bins_array)
a = []
b = []
for i in xrange(bins_array.shape[0]):
    bulged_size = TF_array_bulged[i].shape[0]
    bulgeless_size = TF_array_bulgeless[i].shape[0]
    if np.logical_and(bulged_size>0,bulgeless_size>0):
        if bulgeless_size<bulged_size:
            bulgeless = TF_array_bulgeless[i]
            bulged = np.random.choice(TF_array_bulged[i],
                                      size = (bulgeless_size))
            a.append(bulged)
            b.append(bulgeless)
        else:
            bulged = TF_array_bulged[i]
            bulgeless = np.random.choice(TF_array_bulgeless[i],
                                         size = (bulged_size))
            a.append(bulged)
            b.append(bulgeless)    

n = []   
m = []
o = []   
q = []
for i in xrange(len(b)):
    for j in (np.asarray(b)[i]):
        n.append(R[red_bulgeless][j])
        m.append(uminusr[red_bulgeless][j])
for i in xrange(len(a)):
    for j in (np.asarray(a)[i]):
        o.append(R[red_bulged][j])
        q.append(uminusr[red_bulged][j])
        
redbulgeless = np.asarray(n)
uminusr_bulgeless = np.asarray(m)
redbulged = np.asarray(o)
uminusr_bulged = np.asarray(q)

plt.figure() #figure 4
ax = plt.gca()
ax.invert_xaxis()
plt.xlim(-16,-25)
plt.ylim(0.5,4)
plt.plot(redbulged,uminusr_bulged, "ro", label = "bulged")
plt.plot(redbulgeless,uminusr_bulgeless, "bo",label = "bulgeless")
plt.plot(g,gv,'g-', linewidth = 3, label = "green valley")
plt.legend(loc="best")
plt.xlabel("Absolute R")
plt.ylabel("U - R")
