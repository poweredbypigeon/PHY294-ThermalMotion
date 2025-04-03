import numpy as np 

t = 60 #seconds, 120 frames, 2 frames per second
r_bee = 9.5e-7 #diamater of the bead
temp = 296.5
temp_diff = 296.5 - 293.15 #difference from 20 degrees C to ambient temp
eta = (1 - 0.02*(temp_diff))*0.1 #in pascal seconds, converted from centipoise
gamma = 6*np.pi*eta*r_bee
print(gamma)