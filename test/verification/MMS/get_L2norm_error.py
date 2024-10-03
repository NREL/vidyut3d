import yt
from sys import argv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def get_oned_data(arr3d,axialdir):
    res=np.array(arr3d.shape)
    mid=res/2
    mid=mid.astype(int)
    linedata=np.zeros(res[axialdir])
    trans1dir=(axialdir+1)%3
    trans2dir=(axialdir+2)%3
    
    if(axialdir==0):
        linedata=arr3d[:,mid[trans1dir],mid[trans2dir]]
    if(axialdir==1):
        linedata=arr3d[mid[trans2dir],:,mid[trans1dir]]
    if(axialdir==2):
        linedata=arr3d[mid[trans1dir],mid[trans2dir],:]

    return(linedata)
    

plt.rcParams['font.size'] = 16
ds=yt.load(argv[1])
prob_lo=ds.domain_left_edge.d
prob_hi=ds.domain_right_edge.d

minlev=0
maxlev=ds.index.max_level
lengths=prob_hi-prob_lo
ncells=ds.domain_dimensions
axialdir=np.argmax(lengths)
axialdir_char=chr(ord('x')+axialdir)
covgrid_lev=maxlev
res=np.array([ncells[0]* (2**covgrid_lev),ncells[1]* (2**covgrid_lev),ncells[2]* (2**covgrid_lev)])
dx_frb=lengths/res
fields_load=["AR","Potential"]
ad = ds.covering_grid(level=covgrid_lev, left_edge=prob_lo, dims=res, fields=fields_load)
pot=np.array(ad["Potential"])
s1=np.array(ad["AR"])

pot_1d=get_oned_data(pot,axialdir)
s1_1d=get_oned_data(s1,axialdir)

x=np.linspace(prob_lo[axialdir]+0.5*dx_frb[axialdir],\
        prob_hi[axialdir]-0.5*dx_frb[axialdir],res[axialdir])
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
#fig.suptitle("AR and potential solution along "+dir_char+" direction ")
plt.tight_layout()
plt.savefig("err_"+dir_char+".png")
#=======================================

