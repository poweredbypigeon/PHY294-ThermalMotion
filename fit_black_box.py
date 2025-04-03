# -*- coding: utf-8 -*-
"""
This program will find the best fit of a given function to a given set
of data (including errorbars). It prints the results, with uncertainties.
Then it plots the graph and displays it to the screen, and also saves
a copy to a file in the local directory. Below the main graph is a 
residuals graph, the difference between the data and the best fit line.

There is also a function which will load data from a file. More convenient.
The first line of the file is ignored (assuming it's the name of the variables).
After that the data file needs to be formatted: 
number space number space number space number newline
Do NOT put commas in your data file!! You can use tabs instead of spaces.
The data file should be in the same directory as this python file.
The data should be in the order:
x_data y_data x_uncertainty y_uncertainty
"""


import scipy.optimize as optimize
import numpy as np
import matplotlib.pyplot as plt
from pylab import loadtxt

import inspect

def load_data(filename):
    data=loadtxt(filename, usecols=(0,1,2,3), skiprows=1, unpack=True)
    return data


def plot_fit(my_func, xdata, ydata, xerror=None, yerror=None, init_guess=None, font_size=14, 
             xlabel="Independant Variable (units)", ylabel="Dependent Variable (units)", title="Generic title", log_scale="", filename="graph.pdf"):    
    # xdata and ydata needs to be in np array form
    plt.rcParams.update({'font.size': font_size})
    plt.rcParams['figure.figsize'] = 10, 9
    # Change the fontsize of the graphs to make it easier to read.
    # Also change the picture size, useful for the save-to-file option.
               
    popt, pcov = optimize.curve_fit(my_func, xdata, ydata, sigma=yerror, p0=init_guess, maxfev=100000)
    # The best fit values are popt[], while pcov[] tells us the uncertainties.
    print(tuple(popt))

    puncert = np.sqrt(np.diagonal(pcov))
    # The uncertainties are the square roots of the diagonal of the covariance matrix
    
    print("Best fit parameters, with uncertainties, but not rounded off properly:")
    for i in range(len(popt)):
        print(popt[i], "+/-", puncert[i])
    
    start = min(xdata)
    stop = max(xdata)    
    xs = np.arange(start,stop,(stop-start)/1000) 
    curve = my_func(xs, *popt) 
    # (x,y) = (xs,curve) is the line of best fit for the data in (xdata,ydata).
    # It has 1000 points to make it look smooth.
    # Note: the "*" tells Python to send all the popt values in a readable way.
    
    fig, (ax1,ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 1]})
    # Make 2 graphs above/below each other: ax1 is top, ax2 is bottom.
    # The gridspec_kw argument makes the top plot 2 times taller than the bottom plot.
    # You can adjust the relative heights by, say, changing [2, 1] to [3, 1].
    
    ax1.set_title(title, fontsize=20)
    
    ax1.errorbar(xdata, ydata, yerr=yerror, xerr=xerror, fmt=".", label="data", color="black", ecolor="black", lw=1)
    # normally ecolor is blue
    
    # to check the accuracy of the fit, https://jakevdp.github.io/PythonDataScienceHandbook/04.03-errorbars.html
    # changed the data points to black and best fit graph to red. Since I didn't know if the peaks were error bars I set the error color to blue.
    
    # ax1.errorbar(xdata, ydata, yerr=0, xerr=0, fmt=".", label="data", color="black", ecolor="lightgray", lw=1)
    # yerror, xerror temporarily set to 0
    # Plot the data with error bars, fmt makes it data points not a line, label is
    # a string which will be printed in the legend, you should edit this string.
    
    parameters = inspect.signature(my_func).parameters
    parameter_names = list(parameters.keys())

    plot_label = "best fit: "
    for i in range(1, len(parameter_names)):
        plot_label += parameter_names[i] + ": "
        plot_label += str(round(tuple(popt)[i-1], 3))
        # figure out MatPlotLib popt and what it even is
        plot_label += " "

    ax1.plot(xs, curve, label=plot_label, color="red")
    # make it suitable for all arguments, 5.3f seems to be to 3 floating pts
    # ax1.plot(xs, [my_func(x, init_guess[0], init_guess[1]) for x in xs], label="initial guess", color="blue")
    # Plot the best fit curve on top of the data points as a line.
    # NOTE: you may want to change the value of label to something better!!

    ax1.legend(loc='upper right')
    # Prints a box using what's in the "label" strings in the previous two lines.
    # loc specifies the location
    
    #ax1.axhline(y=4.15759153, color="black")

    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    # label the axes
    
    log_scale = log_scale.lower()
    if log_scale == "x":
        ax1.set_xscale('log')
    elif log_scale == "y":
        ax1.set_yscale('log')
    elif log_scale == "xy" or log_scale == "yx":
        ax1.set_xscale('log')
        ax1.set_yscale('log')
    # uncomment out the above two lines if you want to make it log-log scale
    
    residual = ydata - my_func(xdata, *popt)
    ax2.errorbar(xdata, residual, yerr=yerror, xerr=xerror, fmt=".", color="black", lw=1)
    # Plot the residuals with error bars.
    
    ax2.axhline(y=0, color="black")    
    # Plot the y=0 line for context.
    
    
    
    ax2.set_xlabel(xlabel)
    ax2.set_ylabel("Residuals")
    ax2.set_title("Residuals of the fit")
    # Here is where you change how your graph is labelled.

    fig.tight_layout()
    # Does a decent cropping job of the two figures.
    
    y_mean = np.mean(ydata)
    ss_total = np.sum((ydata - y_mean) ** 2)
    ss_residual = np.sum(residual ** 2)
    r_squared = 1 - (ss_residual / ss_total)
    # https://saturncloud.io/blog/quantifying-the-quality-of-curve-fit-using-python-scipy/
    
    plt.annotate("r^2 = {:.3f}".format(r_squared), (0, 1))
    
    print("R-squared:", r_squared)
    
    plt.show()
    # Show the graph on your screen.

    plt.savefig(filename)  
    # This saves the graph as a file, which will get overwritten
    # every time you run this program, so rename the file if you
    # want to keep multiple files!

    return popt, puncert # return the best fit parameters for further use
    

'''
Todo:
add reduced chi square code somewhere, too. 

'''