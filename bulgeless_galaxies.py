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
def TF_array(data,column,bins_array):
    c = []
    for i in xrange(bins_array.shape[0]):
        #digitize calculates hist bins, in1d identifies where two arrays have 
        #same value
        z = np.where(np.in1d(bins_array[np.digitize(data[...,column],
                                                    bins_array)],
                                                    bins_array[i]))
        c.append(z[0])
    return np.asarray(c)
def index(bin1,galaxy,a,b,c,d,e,f,g,h,k,l,m,r_band,u_band,NUV,NUV_err,u_err,r_err,
          dr8,ra,dec,redshift,fraction):
    for i in xrange(len(bin1)):
        for j in (np.asarray(bin1)[i]):
            a.append(galaxy[j,r_band])
            b.append(galaxy[j,u_band]-galaxy[j,r_band])
            c.append(-1*((-1*galaxy[j,NUV])-galaxy[j,u_band]))
            d.append(np.sqrt((galaxy[j,NUV_err])**2 + (galaxy[j,u_err])**2))
            e.append(np.sqrt((galaxy[j,u_err])**2+(galaxy[j,r_err])**2))
            f.append(galaxy[j,dr8])
            g.append(galaxy[j,ra])
            h.append(galaxy[j,dec])
            k.append(galaxy[j,redshift]) 
            l.append(galaxy[j,u_band])
            m.append(galaxy[j,fraction])
    bulgeless_final = np.column_stack((np.asarray(b), np.asarray(e),
                                      np.asarray(c), np.asarray(d),
                                      np.asarray(k), np.asarray(f),
                                      np.asarray(g), np.asarray(h)))
    red_new = np.asarray(a)
    u_new = np.asarray(l)
    fraction_new = np.asarray(m)
    return fraction_new,u_new,red_new,bulgeless_final
#Data Input####################################################################
read_data = False
if read_data:
    data = np.genfromtxt('combined.csv', delimiter = ',',skip_header = 1)
    np.save('data', data)
    #importing NUV data from bulged/bulgeless galaxies
    NUV_bulged = np.genfromtxt("bulged.csv", delimiter = ",", skiprows = 1)
    np.save('bulged', NUV_bulged)
    NUV_bulgeless = np.genfromtxt("bulgeless.csv", delimiter = ",",
                                  skiprows = 1)
    np.save('bulgeless', NUV_bulgeless)
data = np.load('data.npy')
'''
Column 1 = dr8  ;   Column 3 = ra_1     ;   Column 4 = dec
Column 21 = disk    ;   Column 43 = face on     ;   Column 69 = no bulge
Column 75 = just noticeable     ;   Column 81 = obvious     
Column 85 = dominant    ;   Column 278 = redshift   ;   Column 280 = U
Column 281 = G  ;   Column 282 = R  ;   Column 285 = error in U
Column 287 = error in R
'''
#columns
r_band = 282
r_err = 287
u_band = 280
u_err = 285
g_band = 281
redshift = 278
NUV = 304
NUV_err = 305
dr8 = 1
ra = 3
dec = 4
bulgeless = 69
dominant1 = 85
#selecting data ###############################################################
disk = data[...,21]
face_on = data[...,43]
no_bulge = data[...,bulgeless]
just_noticeable = data[...,75]
obvious = data[...,81]
dominant = data[...,dominant1]
U = data[...,u_band]
G = data[...,g_band]
R = data[...,r_band]
redshift = data[...,redshift]
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
NUV_bulged = np.load('bulged.npy')
NUV_bulgeless =  np.load('bulgeless.npy')
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
TF_array_bulgeless = TF_array(NUV_bulgeless,282,bins_array)  
#extracting bulged data from bins #############################################
TF_array_bulged = TF_array(NUV_bulged,282,bins_array)  
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
#organising indices############################################################
n = []   
m = []
o = []  
p = [] 
q = []
r = []
s = []
t = []
u = []
v = []
w = []
bulgeless,u_bulgeless,redbulgelessfinal, bulgeless_final = index(b,NUV_bulgeless,n,m,o,p,
                                                       q,r,s,t,u,v,w,r_band,
                                                       u_band,NUV,NUV_err,
                                                       u_err,r_err,dr8,ra,dec,
                                                       redshift,no_bulge)      
#f = open('list_commands_bulgeless.txt', 'a')
#b = bulgeless_final
#for n in range(len(b)):
#    f.write('python starpy.py '+str(b[n,0])+' '+str(b[n,1])+' '+str(b[n,2])+' '+str(b[n,3])+' '+str(b[n,4])+'  '+str(b[n,5].astype(int))+' '+str(b[n,6])+' '+str(b[n,7])+'\n')
#f.close()                                                    
n1 = []   
m1 = []
o1 = []  
p1 = [] 
q1 = []
r1 = []
s1 = []
t1 = []
u1 = []       
v1 = []
w1 = []
                                
bulged,u_bulged,redbulgedfinal, bulged_final = index(a,NUV_bulged,n1,m1,o1,p1,q1,r1,
                                              s1,t1,u1,v1,w1,r_band,u_band,NUV,
                                              NUV_err,u_err, r_err,dr8,ra,dec,
                                              redshift,dominant)
#f = open('list_commands_bulged.txt', 'a')
#b = bulged_final
#for n in range(len(b)):
#    f.write('python starpy.py '+str(b[n,0])+' '+str(b[n,1])+' '+str(b[n,2])+' '+str(b[n,3])+' '+str(b[n,4])+'  '+str(b[n,5].astype(int))+' '+str(b[n,6])+' '+str(b[n,7])+'\n')
#f.close()                                               
#test plotting - R,U-R#########################################################
plt.figure() #figure 4 U-R,R
ax = plt.gca()
ax.invert_xaxis()
plt.xlim(-16,-25)
plt.ylim(0.5,4)
plt.plot(redbulgedfinal,bulged_final[...,0], "ro", label = "bulged")
plt.plot(redbulgelessfinal,bulgeless_final[...,0], "bo",label = "bulgeless")
plt.plot(g,gv,'g-', linewidth = 3, label = "green valley")
plt.legend(loc="best")
plt.xlabel("R")
plt.ylabel("U - R")
#test plotting - NUV-U,U-R#####################################################
plt.figure() #figure 5 NUV-U,U-R
ax = plt.gca()
ax.invert_xaxis()
plt.plot(bulged_final[...,0],bulged_final[...,2], "ro", label = "bulged")
plt.plot(bulgeless_final[...,0],bulgeless_final[...,2], "bo",
         label = "bulgeless")
#plt.plot(g,gv,'g-', linewidth = 3, label = "green valley")
plt.legend(loc="best")
plt.ylabel("NUV-U")
plt.xlabel("U - R")
#test plotting - vote fraction vs redshift#####################################
plt.figure() #figure 6 vote fraction, redshift
ax = plt.gca()
plt.plot(redbulgedfinal,bulged_final[...,0], "ro", label = "bulged")
plt.plot(redbulgelessfinal,bulgeless_final[...,0], "bo",label = "bulgeless")
plt.plot(g,gv,'g-', linewidth = 3, label = "green valley")
plt.legend(loc="best")
plt.xlabel("R")
plt.ylabel("U - R")