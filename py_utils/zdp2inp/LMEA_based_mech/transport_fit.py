#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 16:27:42 2023

@author: taareshtaneja
"""

# -*- coding: utf-8 -*-
"""
Created on Mon May  1 16:02:57 2023

@author: taare
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


# This function fit works for both mobility and diffusivity for Dujko et al's N2 streamer case
# Use the exact same of input (x,y) data transformations for both mu and D.

def Arrhenius_fit(x,A,B,C,D,E,F,G,H,I):
    k = A + B*x + C*x**2 + D*x**3 + E*x**4 + F*np.log(x) + G*np.sin(x) + H*np.cos(x) + I*np.exp(-1*x)
    # k = np.exp(A + B*np.log(x) + C/x + D/x**2 + E/x**3)+ F*np.cos(x)
    # k = A + B*x + C*x**2 + D*np.tan(x) + E*x**3 + F*np.cos(G*x) + H*np.sin(I*x) + J*np.log(x)
    return k

data = np.loadtxt("input4.csv",delimiter=",")

xdata = np.log10(data[:,0]+1.01)
ydata = data[:,1] / 1e24


# rate coeffs
parameters, cov = curve_fit(Arrhenius_fit,xdata,ydata)
A,B,C,D,E,F,G,H,I = parameters
print("A,B,C,D,E,F,G = ", parameters)

# # validation
y_validation = Arrhenius_fit(xdata,A,B,C,D,E,F,G,H,I)

# # test
x_min = 0.3
x_max = np.log10(3000)
num_points = 500
# x1 = np.linspace(x_min,x_max,num_points)
x_test = np.linspace(x_min,x_max,num_points)
y_test = Arrhenius_fit(x_test,A,B,C,D,E,F,G,H,I)
# # y_test = (9.0e-09)*(x_test**0.7)*np.exp(-157821.2/x_test) / (eVtoK**0.7)

plt.loglog(10**xdata,1e24*ydata)
plt.loglog(10**xdata,1e24*y_validation)
plt.loglog(10**x_test,1e24*y_test)
# plt.ylim(1e23,1e26)
# 
# # validation error
# plt.plot(xdata,(ydata-y_validation)/ydata)

