import yt
from sys import argv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
    

clength   = 1.0
cwidth    = 0.125
cdepth    = 0.125
ar     = (clength/cwidth)


ds=yt.load(argv[1])
axialdir=np.argmax(ds.domain_dimensions)
axialdir_char=chr(ord('x')+axialdir)
res=100
slicedepth = cdepth/2
slc = ds.slice((axialdir+2)%3,slicedepth)
frb = slc.to_frb((1.0,'cm'),res)
x = np.linspace(0,clength,res)
fld_S1 = np.array(frb["S1"])[res//2,:]
fld_S2 = np.array(frb["S2"])[res//2,:]

c=1.0;
d=0.1;
exactsoln=(np.exp(c*x/d)-1)/(np.exp(c/d)-1)
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

ax[1][0].plot(x,0.5*(x-x**2),'r-',label="Exact solution")
ax[1][0].plot(x,fld_S2,'k*',label="echemAMR",markersize=2)
ax[1][0].legend(loc="best")

im=ax[1][1].imshow(np.array(frb["S2"]),origin="lower")
fig.colorbar(im, ax=ax[1][1])

dir_char=axialdir_char
fig.suptitle("S1 and S2 solution along "+dir_char+" direction ")
plt.savefig("species_"+dir_char+".png")
plt.show()
#=======================================

