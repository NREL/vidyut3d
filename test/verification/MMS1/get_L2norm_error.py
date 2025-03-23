import yt
from sys import argv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['font.size'] = 16
ds=yt.load(argv[1])
axialdir=int(argv[2])
prob_lo=ds.domain_left_edge.d
prob_hi=ds.domain_right_edge.d
maxlev=ds.index.max_level
lengths=prob_hi-prob_lo
mids=0.5*(prob_lo+prob_hi)
ncells=ds.domain_dimensions
axialdir_char=chr(ord('x')+axialdir)
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

pot_1d=lb["Potential"].value[xarr<1.0]
s1_1d=lb["NI"].value[xarr<1.0]
x=xarr[xarr<1.0]

exactsoln_s1=x**2
exactsoln_pot=(x**4-x)/12.0
err_s1=np.sqrt(np.mean((s1_1d-exactsoln_s1)**2))
err_pot=np.sqrt(np.mean((pot_1d-exactsoln_pot)**2))
print(dx_frb[axialdir],err_s1,err_pot)
#=======================================

#=======================================
#Plot solutions
#=======================================
fig,ax=plt.subplots(1,2,figsize=(8,4))
ax[0].plot(x,exactsoln_s1,'r-',label="Exact solution")
ax[0].plot(x,s1_1d,'k*',label="Computed",markersize=2)
ax[0].legend(loc="best")

ax[1].plot(x,(x**4-x)/12.0,'r-',label="Exact solution")
ax[1].plot(x,pot_1d,'k*',label="Computed",markersize=2)
ax[1].legend(loc="best")

dir_char=axialdir_char
#fig.suptitle("NI and potential solution along "+dir_char+" direction ")
plt.tight_layout()
plt.savefig("err_"+dir_char+".png")
#=======================================

