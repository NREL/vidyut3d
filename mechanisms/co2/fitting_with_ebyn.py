import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from sys import argv
import sys
from io import StringIO

def ebyn_fit1(x,A,B,C,D,E):
    f = A*np.log(x) + B + C*x + D*x**2 + E*x**3
    return(f)

def ebyn_fit2(x,A,B,C,D,E):
    f = A*np.log(x) + B + C/x + D/x**2 + E/x**3
    return(f)

def ebyn_fit3(x,A,B,C,D,E):
    f = A*E+B*np.tanh(C*np.log(x)+D)
    return(f)

def ebyn_fit3_string(A,B,C,D,E,scaling,is_exp):
    string = str(A)+"*"+str(E)+"+"+str(B)+"*std::tanh"+"("+str(C)+"*log(ebyn)+"+str(D)+")"
    if(is_exp==True):
        string = "std::exp("+string+")"
    string = str(scaling)+"*"+string
    return(string)

def ebyn_fit2_string(A,B,C,D,E,scaling,is_exp):
    
    string = str(A)+"*std::log(ebyn)+"+str(B)+"+"+str(C)+"/std::log(ebyn)+"+str(D)+"/std::pow(ebyn,2.0)+"+str(E)+"/std::pow(ebyn,3.0)"
    if(is_exp==True):
        string = "std::exp("+string+")"
    string = str(scaling)+"*"+string
    return(string)


def ebyn_fit1_string(A,B,C,D,E,scaling,is_exp):
    string = str(A)+"*std::log(ebyn)+"+ str(B)+"+"+str(C)+"*ebyn+"+ str(D)+"*std::pow(ebyn,2.0)+"+str(E)+"*std::pow(x,3.0)"
    if(is_exp==True):
        string = "std::exp("+string+")"
    string = str(scaling)+"*"+string
    return(string)

def is_number(word):
    try:
        float(word)
    except:
        return(False)
    return(True)

def is_line_of_numbers(string):
    splt=string.split()
    for i in range(len(splt)):
        if(is_number(splt[i])==False):
            return(False)
    return(True)


bolsigoutfile=argv[1]
search=argv[2]
try:
    fit_type=int(argv[3])
except:
    fit_type=3
infile=open(bolsigoutfile,'r')
lines=infile.readlines()
nlines=len(lines)

for l,lineno in  zip(lines,range(nlines)):
    if(search in l):
        break

all_lines=""
for i in range(lineno,nlines):
    line=lines[i]
    splt=line.split()
    if(len(splt)>0):
        if(is_line_of_numbers(line)):
            all_lines=all_lines+line
            #print(np.loadtxt(StringIO(line)))
    else:
        break

arr=np.loadtxt(StringIO(all_lines))
print(arr.shape)

minval=np.min(arr[np.nonzero(arr[:,1]),1])
for i in range(len(arr[:,1])):
    if(arr[i,1]==0.0):
        arr[i,1]=minval

avgval=np.mean(arr[:,1])

ebyn_fit=ebyn_fit1
ebyn_fit_string=ebyn_fit1_string
if(fit_type==1):
    ebyn_fit=ebyn_fit1
    ebyn_fit_string=ebyn_fit1_string
elif(fit_type==2):
    ebyn_fit=ebyn_fit2
    ebyn_fit_string=ebyn_fit2_string
elif(fit_type==3):
    ebyn_fit=ebyn_fit3
    ebyn_fit_string=ebyn_fit3_string
else:
    print("Unknown fit type")
    sys.exit()

parameters, cov = curve_fit(ebyn_fit,arr[:,0],np.log(arr[:,1]/avgval))
A,B,C,D,E=parameters
fitfunc=np.exp(ebyn_fit(arr[:,0],A,B,C,D,E))*avgval
fitstring=ebyn_fit_string(A,B,C,D,E,avgval,True)
print(parameters)
print(fitstring)
np.savetxt("fitting_"+search.replace(" ","")+".dat",np.transpose(np.vstack((np.transpose(arr),fitfunc))))
