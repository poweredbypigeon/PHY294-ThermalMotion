import numpy as np
import matplotlib.pyplot as plt
import os

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

t = 60 #seconds, 120 frames, 2 frames per second
r_bee = 9.5e-7 #diamater of the bead
temp = 296.5
temp_diff = 296.5 - 293.15 #difference from 20 degrees C to ambient temp
eta = (1 - 0.02*(temp_diff))*0.1 #in pascal seconds, converted from centipoise

gamma = 6*np.pi*eta*r_bee

xdata = [] #displacement in x in pixels
ydata = [] #displacement in y in pixels

for i in range(len(filelist)):
    xdata.append(np.loadtxt(filelist[i], skiprows = 2, usecols = 0))
    ydata.append(np.loadtxt(filelist[i], skiprows = 2, usecols = 1))
    
xdiff = np.zeros((len(xdata), 119))
ydiff = np.zeros((len(xdata), 119))*0.1155e-6 #correction factor, pixels to m

for a in range(len(xdata)):
    for b in range(len(xdata[a])-1):
        xdiff[a][b]=xdata[a][b+1]-xdata[a][b] 
        ydiff[a][b]=xdata[a][b+1]-xdata[a][b]

xdiff_corrected = xdiff*0.1155e-6 #correction factor; pixels to m
ydiff_corrected = ydiff*0.1155e-6 #correction factor; pixels to m

meansquarex = np.average(xdiff_corrected**2)
meansquarey = np.average(ydiff_corrected**2)

D = (meansquarex + meansquarey)/(6*t)

k = D*gamma/temp 

plt.plot(t, meansquarex)

step_size = np.sqrt(xdiff_corrected**2 + ydiff_corrected**2)
unique_values = np.unique(step_size)
av_unique_values = average(step_size) #expected average emission
std_unique_values = stdev(step_size)
plt.figure(0)
plt.hist(step_size, bins=unique_values)