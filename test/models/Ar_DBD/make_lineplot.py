import yt
from sys import argv
import matplotlib.pyplot as plt
import numpy as np
import glob
import matplotlib as mpl
from matplotlib.collections import LineCollection

fn_pattern= argv[1]
varname=argv[2]
fn_list=[]
try:
    fn_list = sorted(glob.glob(fn_pattern), key=lambda f: int(f.split("plt")[1]))
except:
    if(fn_list==[]):
        using_file=False
        try:
            print("using file of plotfiles..")
            infile=open(argv[1],'r')
            for line in infile:
                fn_list.append(line.split()[0])
            infile.close()
        except:
            fn_list.append(argv[1])

print(fn_list)
if(len(argv) > 3):
    minval=float(argv[2])
    maxval=float(argv[3])
    set_minmax=True


time_min=yt.load(fn_list[0]).current_time
time_max=yt.load(fn_list[-1]).current_time

ds=yt.load(fn_list[0])
prob_lo=ds.domain_left_edge.d
prob_hi=ds.domain_right_edge.d
probsize=prob_hi-prob_lo
ncells=ds.domain_dimensions
res=ncells[0]
dx=probsize[:]/ncells[:]
xarr = np.linspace(prob_lo[0]+0.5*dx[0],prob_hi[0]-0.5*dx[0],res)

eps=1e-12
fig, ax = plt.subplots(1,1,figsize=(8,5))
lineset=[]
maxfieldval=-1e50
minfieldval=1e50
for i, fn in enumerate(fn_list):
    ds=yt.load(fn)
    time=ds.current_time
    lb = yt.LineBuffer(ds, (prob_lo[0], (0.5+eps)*(prob_lo[1]+prob_hi[1]), (0.5+eps)*(prob_lo[2]+prob_hi[2])), \
        (prob_hi[0], (0.5+eps)*(prob_lo[1]+prob_hi[1]), (0.5+eps)*(prob_lo[2]+prob_hi[2])), res)
    fieldarr=lb[varname]
    lineset.append(np.transpose(np.vstack((xarr,fieldarr))))
    if(np.max(fieldarr)>maxfieldval):
        maxfieldval=np.max(fieldarr)
    if(np.min(fieldarr)<minfieldval):
        minfieldval=np.min(fieldarr)


ax.set_xlim(prob_lo[0],prob_hi[0])
ax.set_ylim(minfieldval,maxfieldval)
line_segments = LineCollection(lineset,cmap='coolwarm',linestyles='solid')
line_segments.set_array(np.linspace(time_min,time_max,len(fn_list)))
ax.add_collection(line_segments)
axcb = fig.colorbar(line_segments)
axcb.set_label('time (sec)')
plt.sci(line_segments)  # This allows interactive changing of the colormap.
plt.tight_layout()
plt.savefig("lsetplot.png")
plt.show()
