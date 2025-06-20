import yt
from sys import argv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
    

ds=yt.load(argv[1])
prob_lo=ds.domain_left_edge.d
prob_hi=ds.domain_right_edge.d
probsize=prob_hi-prob_lo
ncells=ds.domain_dimensions
mids=0.5*(prob_lo+prob_hi)

axialdir=np.argmax(ds.domain_dimensions)
axialdir_char=chr(ord('x')+axialdir)

minlev=0
maxlev=ds.index.max_level
lengths=prob_hi-prob_lo
res=np.array([ncells[0]* (2**maxlev),ncells[1]* (2**maxlev),ncells[2]* (2**maxlev)])
dx_frb=lengths/res
low_end=[mids[0],mids[1],mids[2]]
high_end=[mids[0],mids[1],mids[2]]
low_end[axialdir]=prob_lo[axialdir]+0.5*dx_frb[axialdir]
high_end[axialdir]=prob_hi[axialdir]-0.5*dx_frb[axialdir]
lb = yt.LineBuffer(ds, tuple(low_end), \
    tuple(high_end), res[axialdir])

xarr=np.linspace(prob_lo[axialdir]+0.5*dx_frb[axialdir],\
            prob_hi[axialdir]-0.5*dx_frb[axialdir],res[axialdir])
fld_S1 = lb["S1"].value[xarr>0.0]
fld_S2 = lb["S2"].value[xarr>0.0]
x=xarr[xarr>0.0]

c=1.0;
d=0.1;
exactsoln=(np.exp(c*x/d)-1)/(np.exp(c/d)-1)
#=======================================

#=======================================
#Plot solutions
#=======================================
fig,ax=plt.subplots(1,2,figsize=(8,4))
ax[0].plot(x,exactsoln,'r-',label="Exact solution")
ax[0].plot(x,fld_S1,'k*',label="echemAMR",markersize=2)
ax[0].legend(loc="best")

ax[1].plot(x,0.5*(x-x**2),'r-',label="Exact solution")
ax[1].plot(x,fld_S2,'k*',label="echemAMR",markersize=2)
ax[1].legend(loc="best")

dir_char=axialdir_char
fig.suptitle("S1 and S2 solution along "+dir_char+" direction ")
plt.savefig("species_"+dir_char+".png")
plt.show()
#=======================================

