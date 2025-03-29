# -*- coding: utf-8 -*-
"""
General Code and Given Constants for the Lab
"""
import numpy as np
import scipy.optimize as opt
import scipy.stats as sp
import matplotlib.pyplot as plt
import os

#propagation of error function
def adsubprop(err_x, err_y): #err_x is the error in the argument x
    err_adsub=np.sqrt(err_x**2+err_y**2)
    return(err_adsub)
    
def muldivprop(x, y, z, err_x, err_y):
    err_muldiv=z*(np.sqrt((err_x/x)**2+(err_y/y)**2))
    return(err_muldiv)

def exprop(exponent, x, err_x):
    err_exp=exponent*(x**exponent)*(err_x/x)
    return(err_exp)
    
t = 60 #seconds, 120 frames, 2 frames per second
r_bee = 9.5e-7 #diamater of the bead
temp = 296.5
temp_diff = 296.5 - 293.15 #difference from 20 degrees C to ambient temp
eta = (1 - 0.02*(temp_diff))*0.1 #in pascal seconds, converted from centipoise
gamma = 6*np.pi*eta*r_bee
"""
Code to Clean and Compile Data
"""
#importing the data collected from the observations of motion.

filelist = []

i = 0
longest = 0
#a placeholder array for data, shape of the number of entries we want by the number of measurements taken

#makes a list of the paths of each file in the desired directory
for root, dirs, files in os.walk("Capture Data"):
    for file in files:
        #takes only files of type txt
        if file.endswith(".txt"):
            i=i+1
            filelist.append(os.path.join(root, file))

xdata = [] #displacement in x in pixels
ydata = [] #displacement in y in pixels

for i in range(len(filelist)):
    xdata.append(np.loadtxt(filelist[i], skiprows = 2, usecols = 0))
    ydata.append(np.loadtxt(filelist[i], skiprows = 2, usecols = 1))
    
xdiffs = np.zeros((len(xdata), 119))
ydiffs = np.zeros((len(xdata), 119))*0.1155e-6 #correction factor, pixels to m

for a in range(len(xdata)):
    for b in range(len(xdata[a])-1):
        xdiffs[a][b]=xdata[a][b+1]-xdata[a][b] 
        ydiffs[a][b]=xdata[a][b+1]-xdata[a][b]

xdiff_reshaped = xdiffs.reshape(1, 2618) #2618 = 22*119
ydiff_reshaped = ydiffs.reshape(1, 2618) #done to flatten array

xdiff = np.round(xdiff_reshaped, 1)
ydiff = np.round(ydiff_reshaped, 1)

xdiff_corrected = xdiff*0.1155e-6 #correction factor; pixels to m
ydiff_corrected = ydiff*0.1155e-6 #correction factor; pixels to m


"""
Code for Propagation of Error
"""
xdata_er = np.ones(2618)*np.std(xdata)
ydata_er = np.ones(2618)*np.std(ydata)

err_radius = 1e-10 #0.1 um
err_visc = 0.005 #g/m*s
err_temp = 0.5 #k
err_eta = (1 - 0.02*(0.05))*0.1
err_avx = np.std(xdata)
err_avy = np.std(ydata)

#error in the summ of the average, no error from division by scalar
for c in range(2618):
    err_avx = adsubprop(err_avx, xdata_er[c])
    err_avy = adsubprop(err_avy, ydata_er[c])
    
err_disp = adsubprop(err_avx, err_avy)

err_drag = muldivprop(eta, r_bee, err_eta, err_radius, eta*r_bee)



"""
Code for Histogram, Rayleigh, and Maximum Likelihood Calculations
"""

step_size = np.sqrt(xdiff[0]**2 + ydiff[0]**2)
#step_size = np.round(step_size, decimals = 2)

def rayleigh(r, D, t):
    return (r/(2*D*t))*np.exp(-(r**2)/(4*D*t))

def maximum_likelihood(r):
    return np.sum(r**2)/(2*len(r))

#identifying the unique number of incidences of radiation in the time interval
unique_values = np.unique(step_size)
av_unique_values = np.average(step_size) #expected average emission
std_unique_values = np.std(step_size)

#defining the normal distribution
normal = sp.norm.pdf(step_size, loc=av_unique_values, scale=std_unique_values)
p_opt, p_cov = opt.curve_fit(rayleigh, step_size, normal)

rayleigh_k = p_opt[0]*gamma/temp

#Histogram of incidence of counts per time interval
plt.figure(0)
plt.hist(step_size, bins=unique_values)
plt.plot(step_size, normal, label = '')
plt.plot(step_size, (step_size/(maximum_likelihood(step_size)))*np.exp(-(step_size**2)/(2*maximum_likelihood(step_size))), label = '')
plt.xlabel('')
plt.ylabel('')
plt.title('')
plt.figtext(0.05, 0.93, "Figure 1")
plt.legend()
plt.grid()
plt.savefig('')
plt.ylim(0,10)

"""
Code for Einstein Relationship Calculations
"""
meansquarex = np.average(xdiff_corrected**2)
meansquarey = np.average(ydiff_corrected**2)

D = (meansquarex + meansquarey)/(6*t)

einstein_k = D*gamma/temp 