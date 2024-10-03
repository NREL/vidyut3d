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
fintime=float(argv[2])
vel=float(argv[3])
prob_lo=ds.domain_left_edge.d
prob_hi=ds.domain_right_edge.d

minlev=0
maxlev=ds.index.max_level
lengths=prob_hi-prob_lo
midx=0.5*(prob_hi[0]+prob_lo[0])
ncells=ds.domain_dimensions
axialdir=np.argmax(lengths)
axialdir_char=chr(ord('x')+axialdir)
covgrid_lev=maxlev
res=np.array([ncells[0]* (2**covgrid_lev),ncells[1]* (2**covgrid_lev),ncells[2]* (2**covgrid_lev)])
dx_frb=lengths/res
fields_load=["AR"]
ad = ds.covering_grid(level=covgrid_lev, left_edge=prob_lo, dims=res, fields=fields_load)
s1=np.array(ad["AR"])

s1_1d=get_oned_data(s1,axialdir)

x=np.linspace(prob_lo[axialdir]+0.5*dx_frb[axialdir],\
        prob_hi[axialdir]-0.5*dx_frb[axialdir],res[axialdir])

exactsoln_s1=np.sin(np.pi*(x-vel*fintime)-np.sin(np.pi*(x-vel*fintime))/np.pi)
err_s1=np.sqrt(np.mean((s1_1d-exactsoln_s1)**2))
print(dx_frb[axialdir],err_s1)
#=======================================

#=======================================
#Plot solutions
#=======================================
fig,ax=plt.subplots(1,1,figsize=(4,4))
ax.plot(x,exactsoln_s1,'r-',label="Exact solution")
ax.plot(x,s1_1d,'k*',label="Computed",markersize=2)
ax.legend(loc="best")

dir_char=axialdir_char
#fig.suptitle("AR and potential solution along "+dir_char+" direction ")
plt.tight_layout()
plt.savefig("solncompare_"+dir_char+".png")
#=======================================

