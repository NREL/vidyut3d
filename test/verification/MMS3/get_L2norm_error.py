import yt
from sys import argv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

n0=1e6
alpha=1.60217662e-19/8.854187817e-12
ME = 9.10938188e-31
M_AMU = 1.66054e-27
nu=10000.0

plt.rcParams['font.size'] = 16
ds=yt.load(argv[1])
axialdir=int(argv[2])
finaltime=float(argv[3])
dt=float(argv[4])

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

pot_1d=lb["Potential"].value
s1_1d=lb["E"].value
s2_1d=lb["HEp"].value
s3_1d=lb["Electron_energy"].value
x=np.linspace(prob_lo[axialdir]+0.5*dx_frb[axialdir],\
        prob_hi[axialdir]-0.5*dx_frb[axialdir],res[axialdir])
exactsoln_s1=n0*(1.0+np.sin(np.pi*x/lengths[axialdir]))*np.exp(-5.0*finaltime)
exactsoln_s2=n0*(1.0+np.sin(np.pi*x/lengths[axialdir]))*np.exp(-5.0*finaltime)
#lamda_ee=nu*(2.0*ME)/(4.0*M_AMU)
lamda_ee=10.0
exactsoln_s3=n0*(1.0+np.sin(np.pi*x/lengths[axialdir]))*np.exp(-lamda_ee*finaltime)
exactsoln_pot=np.zeros(res[axialdir])
err_s1=np.sqrt(np.mean((s1_1d-exactsoln_s1)**2))
err_s2=np.sqrt(np.mean((s2_1d-exactsoln_s2)**2))
err_s3=np.sqrt(np.mean((s3_1d-exactsoln_s3)**2))
err_pot=np.sqrt(np.mean((pot_1d-exactsoln_pot)**2))
print(dt,err_s1,err_s2,err_s3,err_pot)
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

