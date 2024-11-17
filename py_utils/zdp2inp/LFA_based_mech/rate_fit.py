# -*- coding: utf-8 -*-
"""
Created on Mon May  1 16:02:57 2023

@author: taare
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def Arrhenius_fit(x,A,B,C,D,E,F,G,H,I):
    k = np.exp(A + B*np.log(x) + C/x + D/x**2 + E/x**3)
    # k = A + B*x + C*x**2 + D*x**3 + E*x**4 + F*x**5
    return k

def Arrhenius_fit2(x,A,B,C,D,E,F,G,H,I):
    k = np.exp(A + B*np.log(x) + C/x)
    return k

data = np.loadtxt("input4.csv",delimiter=",")

xdata = data[:,0]
ydata = data[:,1]

# rate coeffs
parameters, cov = curve_fit(Arrhenius_fit2,xdata,ydata)
A,B,C,D,E,F,G,H,I = parameters
print("A,B,C,D,E,F,G = ", parameters)

# # validation
y_validation = Arrhenius_fit2(xdata,A,B,C,D,E,F,G,H,I)

# # test
x_min = 0.01
x_max = xdata[-1]
num_points = 5000
# x1 = np.linspace(x_min,x_max,num_points)
x_test = np.linspace(x_min,x_max,num_points)
y_test = Arrhenius_fit2(x_test,A,B,C,D,E,F,G,H,I)
# # y_test = (9.0e-09)*(x_test**0.7)*np.exp(-157821.2/x_test) / (eVtoK**0.7)

plt.semilogy(xdata,ydata)
# plt.semilogy(xdata,y_validation)
plt.semilogy(x_test,y_test)
plt.ylim(1e-50,1e-5)
plt.xlim(-100,50)
# 
# # validation error
# plt.plot(10**xdata,(10**ydata-10**y_validation)/10**y_validation)

# NH3 E1
# -2.73349619e+01 -2.71647193e-01 -7.81463651e+04 -5.30850417e+08
#   1.00000000e+00
  
# NH3 E2
# -2.78603532e+01 -3.60698526e-02 -1.96732080e+05 -1.28507578e+08
#   1.00000000e+00  

# NH3 ionization
# -3.39902168e+01  3.29753098e-01 -2.06638805e+05 -5.98397173e+08
#   1.00000000e+00

  
  