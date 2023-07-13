import yt
from sys import argv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
    

axialdir = int(argv[2])
clength   = 1.0
cwidth    = 0.125
cdepth    = 0.125
ar     = (clength/cwidth)

axialdir_char=chr(ord('x')+axialdir)

ds=yt.load(argv[1])
res=100
slicedepth = cdepth/2
slc = ds.slice((axialdir+2)%3,slicedepth)
frb = slc.to_frb((1.0,'cm'),res)
x = np.linspace(0,clength,res)
fld_S1 = np.array(frb["S1"])[res//2,:]
fld_pot = np.array(frb["Potential"])[res//2,:]

exactsoln=x**2
#=======================================

#=======================================
#Plot solutions
#=======================================
fig,ax=plt.subplots(2,2,figsize=(8,4))
ax[0][0].plot(x,exactsoln,'r-',label="Exact solution")
ax[0][0].plot(x,fld_S1,'k*',label="echemAMR",markersize=2)
ax[0][0].legend(loc="best")

im=ax[0][1].imshow(np.array(frb["S1"]),origin="lower")
fig.colorbar(im, ax=ax[0][1])

ax[1][0].plot(x,(x**4-x)/12.0,'r-',label="Exact solution")
ax[1][0].plot(x,fld_pot,'k*',label="echemAMR",markersize=2)
ax[1][0].legend(loc="best")

im=ax[1][1].imshow(np.array(frb["Potential"]),origin="lower")
fig.colorbar(im, ax=ax[1][1])

dir_char=chr(ord('x')+int(argv[2]))
fig.suptitle("S1 and S2 solution along "+dir_char+" direction ")
plt.savefig("species_"+dir_char+".png")
plt.show()
#=======================================
