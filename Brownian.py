# -*- coding: utf-8 -*-
"""
General Code and Given Constants for the Lab
"""
import numpy as np
import scipy.optimize as opt
import scipy.stats as sp
import matplotlib.pyplot as plt
import os

#defining the average and standard deviation functions to be used later
def average(container):
    avg = np.sum(container)/len(container)
    return(avg)

def stdev(container):
    mean_x = average(container)
    #an array of the differences between values of the array and the mean of the array
    avdiff_x = container - mean_x
    stdev_x = np.sqrt((sum(avdiff_x**2))/(len(container)-1))
    #tells the function to return the standard deviation of the array
    return(stdev_x)

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

#filelist = []
#
#i = 0
#longest = 0
##a placeholder array for data, shape of the number of entries we want by the number of measurements taken
#
##makes a list of the paths of each file in the desired directory
#for root, dirs, files in os.walk("Capture Data"):
#    for file in files:
#        #takes only files of type txt
#        if file.endswith(".txt"):
#            i=i+1
#            filelist.append(os.path.join(root, file))


SOURCE  = "C:\\Users\\compl\\Documents\\Engineering Science\\Y2S2\\PHY294H1\\Labs\\L3+4\\good_data"

filenames = os.listdir(SOURCE)

filelist = []

for file in filenames:
    filelist.append(os.path.join(SOURCE, file)) 

'''
filelist = np.ndarray([])
for file in filenames:# Assuming to_combine is in current directory.
    for value in np.loadtxt("good_data\\" + file, skiprows=2, unpack=True):
        # print(value)
        filelist = np.append(filelist, value)
'''
# print(filelist)

xdata = [] #displacement in x in pixels
ydata = [] #displacement in y in pixels

for i in range(len(filelist)):
    # print(np.loadtxt(filelist[i], skiprows = 2, usecols = 0))
    xdata.append(np.loadtxt(filelist[i], skiprows = 2, usecols = 0))
    ydata.append(np.loadtxt(filelist[i], skiprows = 2, usecols = 1))



xdiff = np.zeros((len(xdata), 119))
ydiff = np.zeros((len(xdata), 119))*0.1155e-6 #correction factor, pixels to m

for a in range(len(xdata)):
    for b in range(len(xdata[a])-1):
        # print(a, b)
        xdiff[a][b]=xdata[a][b+1]-xdata[a][b] 
        ydiff[a][b]=xdata[a][b+1]-xdata[a][b]

xdiff_corrected = xdiff*0.1155e-6 #correction factor; pixels to m
ydiff_corrected = ydiff*0.1155e-6 #correction factor; pixels to m

step_size = np.sqrt(xdiff_corrected**2 + ydiff_corrected**2)

"""
Code for Propagation of Error
"""
xdata_er = stdev(xdata)
ydata_er = stdev(ydata)
xdiff_er = adsubprop(xdata_er, xdata_er)
ydiff_er = adsubprop(ydata_er, ydata_er)
xdiff_corrected_er = xdiff_er*0.1155e-6
ydiff_corrected_er = ydiff_er*0.1155e-6
#step_size_er = exprop(1/2, xdiff_corrected**2 + ydiff_corrected**2, adsubprop(exprop(2, xdiff_corrected, xdiff_corrected_er), exprop(2, ydiff_corrected, ydiff_corrected_er)))

eta_er = (1 - 0.02*(0.05))*0.1
gamma_er = 6*np.pi*muldivprop(eta, r_bee, eta*r_bee, eta_er, 0.1e-7/2)


"""
Code for Histogram, Rayleigh, and Maximum Likelihood Calculations
"""
def rayleigh(r, D, t):
    return (r/(2*D*t))*np.exp(-(r**2)/(4*D*t))

def maximum_likelihood(r):
    return np.sum(r**2)/(2*len(r))

#identifying the unique number of incidences of radiation in the time interval
unique_values = np.unique(step_size)
av_unique_values = average(step_size) #expected average emission
std_unique_values = stdev(step_size)

#defining the normal distribution
normal = sp.norm.pdf(step_size, loc=av_unique_values, scale=std_unique_values)

# doesn't work
# p_opt, p_cov = opt.curve_fit(rayleigh, step_size, normal)

#rayleigh_k = p_opt[0]*gamma/temp

#Histogram of incidence of counts per time interval
plt.figure(0)
plt.hist(step_size)#, bins=unique_values)
#plt.plot(step_size, normal, label = '')
#plt.plot(step_size, (step_size/(maximum_likelihood(step_size)))*np.exp(-(step_size**2)/(2*maximum_likelihood(step_size))), label = '')
#plt.xlabel('')
#plt.ylabel('')
#plt.title('')
#plt.figtext(0.05, 0.93, "Figure 1")
#plt.legend()
plt.grid()
#plt.savefig('')

"""
Code for Einstein Relationship Calculations
"""
meansquarex = np.average(xdiff_corrected**2)
meansquarey = np.average(ydiff_corrected**2)

D = (meansquarex + meansquarey)/(6*t)

einstein_k = D*gamma/temp 