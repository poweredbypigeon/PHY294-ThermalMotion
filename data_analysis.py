import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd 
import os 
import uncertainties 
from uncertainties import umath
import fit_black_box as bb 
import math
from scipy.optimize import curve_fit


# Constants
r = uncertainties.ufloat(9.5e-7, 5e-8) # bead radius, divided by 2. 
eta = uncertainties.ufloat(1.00, 0.05) # viscosity
T = uncertainties.ufloat(296.5, 0.5) # temperature
eta_adjusted = eta * umath.pow(0.98, T - 293.15) # adjusted viscosity
print("Viscosity adjusted for temperature", eta_adjusted)
gamma = 6 * math.pi * eta_adjusted * r # drag coefficient

print("Calculated Stokes Radius", gamma)

resolution = uncertainties.ufloat(0.1204e-6, 0.003e-6) # camera resolution
# bead = uncertainties.ufloat(2.0e-6, 0.1e-6) # bead size
step_size = uncertainties.ufloat(0.5, 0.03) # step size


LIST_SIZE = 120 # 120 data points

SOURCE  = "C:\\Users\\compl\\Documents\\Engineering Science\\Y2S2\\PHY294H1\\Labs\\L3+4\\good_data"


def adsubprop(series1, series2, unc1, unc2):
    if series1 is None:
        return series2, unc2
    if series2 is None:
        return series1, unc1
    result = series1 + series2  # or series1 - series2
    uncertainty = np.sqrt(unc1**2 + unc2**2)
    return result, uncertainty

# Example: Multiplying two series with uncertainties
def muldivprop(series1, series2, unc1, unc2):
    result = series1 * series2  # or series1 / series2
    relative_unc = np.sqrt((unc1 / series1)**2 + (unc2 / series2)**2)
    uncertainty = result * relative_unc
    return result, uncertainty

# Example: Exponentiation with uncertainty
def exprop(series, exponent, unc):
    result = series**exponent
    relative_unc = np.abs(exponent) * (unc / series)
    uncertainty = result * relative_unc
    return result, uncertainty

def normalize(df):
    
    x_series = df.iloc[:, 0]
    x_series -= x_series[0]
    y_series = df.iloc[:, 1]
    y_series -= y_series[0]
    df = pd.DataFrame({"x": x_series, "y": y_series})
    return df

def conversion(df):
    x_series = df.iloc[:, 0]
    x_series = x_series * resolution.n
    y_series = df.iloc[:, 1]
    y_series = y_series * resolution.n 
    x_unc = pd.Series([1e-7] * x_series.size) # uncertainty in position should be a lot bigger than camera uncertainty
    y_unc = pd.Series([1e-7] * y_series.size)
    df = pd.DataFrame({"x": x_series, "y": y_series, "x_unc": x_unc, "y_unc": y_unc})
    return df

def get_x2_series(df):
    x_series = df.iloc[:, 0]
    x_series -= x_series[0]
    x_series = x_series * x_series 

    return x_series

def linear(x, m, b):
    return m*x + b 


def rayleigh(r, D):
    t = step_size.n
    return (r/(2*D*t))*np.exp(-(r**2)/(4*D*t))

def maximum_likelihood(r):
    return np.sum(r**2)/(2*len(r))

def get_distances(df):
    distances = []
    x_series = df.iloc[:, 0]
    y_series = df.iloc[:, 1]
    for i in range(1, len(x_series)):
        distances.append(math.sqrt((x_series[i] - x_series[i-1])**2 + (y_series[i] - y_series[i-1])**2))
    return distances


filenames = os.listdir(SOURCE)
# everything in good-data folder

# Load the data, skipping the header row

x2_avg = [] # average x^2 values

distances = [] # distances between the two points

s = None 
s_total = None 
s_unc = None


for filename in filenames:
    df = pd.read_csv(os.path.join(SOURCE, filename), delimiter="\t", skiprows=1)
    # print("Number of rows", len(df), "for", filename)
    df = normalize(df)
    df = conversion(df)
    # calculate x^2 and uncertainty 
    x2 = df.iloc[:, 0] * df.iloc[:, 0]
    x2_uncertainty = 2 * abs(df.iloc[:, 0]) * df.iloc[:, 2] # so much bigger than camera uncertainty ignore it
    # x2 = (df.iloc[:, 0] * resolution.n) * (df.iloc[:, 0] * resolution.n)
    # x2_uncertainty = 2 * resolution.n * df.iloc[:, 0] * 1e-7 # so much bigger than camera uncertainty ignore it 

    s_total, s_unc = adsubprop(s_total, x2, s_unc, x2_uncertainty)

    # s = get_x2_series(df) / len(filenames)
    distances.extend(get_distances(df)) # this is just - add the distances to the list
    ''''
    if s_total is None:
        s_total = s
    else:
        s_total += s'
    '''
LIST_SIZE = len(s_total)

s_total = s_total / len(filenames)
s_unc = s_unc / len(filenames)

s_unc[s_unc == 0] = 1e-14 

# print(s_total)
# print(s_unc)
# print(distances)
""" 

guesses = [(4e-13, 0)]

titles = ["X-Squared vs, Time for Beads"]

# for x in [10,20,30,40,50,60,70,80,90,100]:
#     titles.append("Radians as a function of time for an initial angle of 20 degrees and length of " + str(x) + " cm")
#     if x != 70:
#         titles.append("Radians as a function of time for an initial angle of 20 degrees and length of " + str(x) + " cm")

# here this is the theta vs time graph, x = time, y = theta/amplitude

labels = [("Time (s)", "X-Squared (m^2)")]
functions = [linear]

time = np.array([i * step_size.n for i in range(0, LIST_SIZE)])
time_unc = np.array([step_size.s] * LIST_SIZE)

# print(time)
# print(time_unc)

i = 0
init_guess = guesses[i]
font_size = 12
xlabel = labels[i][0]
ylabel = labels[i][1]
title = titles[i]

popt = bb.plot_fit(functions[i], time, s_total, time_unc, s_unc, init_guess=init_guess, font_size=font_size,
            xlabel=xlabel, ylabel=ylabel, title=title, filename=filename.replace(".txt", ".pdf"))
 """
'''
1.0207769625639054e-12 +/- 1.2772905432650444e-14
-8.08275213086696e-13 +/- 1.4848420162119995e-13
R-squared: 0.9627819273532291
'''

# DISTANCES HAS 119 DATA POINTS x 50 FILES = 5950 VALUES 

counts, bin_edges = np.histogram(distances, bins=100, density=True)  # `density=True` normalizes the histogram
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2  # Compute bin centers
print("Bin centers:", bin_centers)
print("Counts:", counts)
# why are the counts so large? they go up to 1M and have fractional values. there are 5950 values 
# EDIT: they're like that as density=True adjusts them so that the integral is 1

# Fit the Gaussian model to the histogram data
popt, pcov = curve_fit(rayleigh, bin_centers, counts, p0=[2.45e-16])  # Initial guess for D, based on calculated viscosity, known T, and k = 1.38 * 10^-23
puncert = np.sqrt(np.diagonal(pcov))
print(f"Fitted parameters")
print(popt)
print(puncert)

# says D = 4.968 * 10^-17. Not sure what pcov is, we need uncertainty values! 

'''
plt.hist(distances, bins=100, color='blue', edgecolor='black')  # Set bins to 20
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Histogram Example')
plt.show()
'''

""" 
def pre_process(df):
    x_series = df.iloc[:, 0]
    x_series -= x_series[0]
    x_series = x_series * x_series 

    y_series = df.iloc[:, 1]
    y_series -= y_series[0]
    y_series = y_series * y_series

    return x_series.tolist(), y_series.tolist()

# Extract columns into lists
# xs = x_series.tolist()
# ys = y_series.tolist()




# print(xs[:5], ys[:5])  # Print first few values to check

Constants-

assume it's a sphere: 
bead radius (r) - 9.5 * 10^-7 +/- 5 * 10^-8 m
viscosity- 1.00 +/- 0.05 g/(cm*s) at 293.15 K, decreases 2% with each degree increase in temp. 
temperature- 296.5 +/- 0.5 K

camera's resolution - 0.1204 +/- 0.003 micrometer/pixel, note that the distances and s_total are in pixels. 
beads - 2.0 +/- 0.1 micrometers (assumed)
step size - 0.5 +/- 0.03 s 
 """