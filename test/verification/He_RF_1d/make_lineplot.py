import yt
from sys import argv
import matplotlib.pyplot as plt
import numpy as np
import glob


fn_pattern= argv[1]
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


for i, fn in enumerate(fn_list):

    ds=yt.load(fn)
    prob_lo=ds.domain_left_edge.d
    prob_hi=ds.domain_right_edge.d
    probsize=prob_hi-prob_lo
    ncells=ds.domain_dimensions
    mids=0.5*(prob_lo+prob_hi)
    
    minlev=0
    maxlev=ds.index.max_level
    lengths=prob_hi-prob_lo
    axialdir=0
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
    xvalid=xarr*(xarr-0.067)
    x=xarr[xvalid<0.0]

    potl=lb["Potential"].value[xvalid<0.0]
    efieldl=lb["Efieldx"].value[xvalid<0.0]
    edenl=lb["E"].value[xvalid<0.0]
    etempl=lb["Electron_Temp"].value[xvalid<0.0]
    eenrgl=lb["Electron_energy"].value[xvalid<0.0]
    iondenl=lb["HEp"].value[xvalid<0.0]
    metadenl=lb["HEm"].value[xvalid<0.0]
    ejl=lb["Electron_Jheat"].value[xvalid<0.0]
    inelhl=lb["Electron_inelasticHeat"].value[xvalid<0.0]
    elhl=lb["Electron_elasticHeat"].value[xvalid<0.0]

    np.savetxt("linedata_"+axialdir_char+"%4.4d.dat"%(i),\
            np.transpose(np.vstack((x,potl,efieldl,edenl,iondenl,metadenl,etempl,eenrgl,ejl,elhl,inelhl))),delimiter="  ")
