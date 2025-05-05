import yt
from sys import argv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

n0=1e6
alpha=1.60217662e-19/8.854187817e-12

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

xvalid=xarr*(xarr-1.0)
x=xarr[xvalid<0.0]
pot_1d=lb["Potential"].value[xvalid<0.0]
s1_1d=lb["E"].value[xvalid<0.0]
s2_1d=lb["HEp"].value[xvalid<0.0]
s3_1d=lb["Electron_energy"].value[xvalid<0.0]
exactsoln_s1=x**3/alpha+n0
exactsoln_s2=x**3/(2.0*alpha)+n0
exactsoln_s3=x**3/(alpha)+n0
exactsoln_pot=(x**5-x)/40.0
err_s1=np.sqrt(np.mean((s1_1d-exactsoln_s1)**2))
err_s2=np.sqrt(np.mean((s2_1d-exactsoln_s2)**2))
err_s3=np.sqrt(np.mean((s3_1d-exactsoln_s3)**2))
err_pot=np.sqrt(np.mean((pot_1d-exactsoln_pot)**2))
print(dx_frb[axialdir],err_s1,err_s2,err_s3,err_pot)
#=======================================

#=======================================
#Plot solutions
#=======================================
fig,ax=plt.subplots(1,3,figsize=(12,4))
ax[0].plot(x,exactsoln_s1,'b-',label="Exact solution (E)")
ax[0].plot(x,s1_1d,'g*',label="Computed (E)",markersize=2)
ax[0].legend(loc="best")

ax[0].plot(x,exactsoln_s2,'r-',label="Exact solution (I)")
ax[0].plot(x,s2_1d,'k*',label="Computed (I)",markersize=2)
ax[0].legend(loc="best")

ax[1].plot(x,exactsoln_pot,'r-',label="Exact solution (Pot)")
ax[1].plot(x,pot_1d,'k*',label="Computed (Pot)",markersize=2)
ax[1].legend(loc="best")

ax[2].plot(x,exactsoln_s3,'r-',label="Exact solution (EEN)")
ax[2].plot(x,s3_1d,'k*',label="Computed (EEN)",markersize=2)
ax[2].legend(loc="best")

dir_char=axialdir_char
plt.tight_layout()
plt.savefig("err_"+dir_char+".png")
#=======================================

