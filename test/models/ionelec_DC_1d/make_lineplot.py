import yt
from sys import argv
import matplotlib.pyplot as plt
import numpy as np
import glob

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
    
    minlev=0
    maxlev=ds.index.max_level
    lengths=prob_hi-prob_lo
    axialdir=np.argmax(lengths)
    axialdir_char=chr(ord('x')+axialdir)

    
    covgrid_lev=maxlev
    res=np.array([ncells[0]* (2**covgrid_lev),ncells[1]* (2**covgrid_lev),ncells[2]* (2**covgrid_lev)])
    dx_frb=probsize/res
    fields_load=["Potential","E","Electron_energy","Electron_Temp","HEp"]
    ad = ds.covering_grid(level=covgrid_lev, left_edge=prob_lo, dims=res, fields=fields_load)

    pot=np.array(ad["Potential"])
    elecden=np.array(ad["E"])
    etemp=np.array(ad["Electron_Temp"])
    eenrg=np.array(ad["Electron_energy"])
    ionden=np.array(ad["HEp"])
    ejheat=np.array(ad["Electron_Jheat"])
    inelheat=np.array(ad["Electron_inelasticHeat"])
    elheat=np.array(ad["Electron_elasticHeat"])
    xarr=np.linspace(prob_lo[axialdir]+0.5*dx_frb[axialdir],\
            prob_hi[axialdir]-0.5*dx_frb[axialdir],res[axialdir])

    potl=get_oned_data(pot,axialdir)
    edenl=get_oned_data(elecden,axialdir)
    etempl=get_oned_data(etemp,axialdir)
    eenrgl=get_oned_data(eenrg,axialdir)
    iondenl=get_oned_data(ionden,axialdir)
    ejl=get_oned_data(ejheat,axialdir)
    elhl=get_oned_data(elheat,axialdir)
    inelhl=get_oned_data(inelheat,axialdir)

    np.savetxt("linedata_"+axialdir_char+"%d.dat"%(i),\
            np.transpose(np.vstack((xarr,potl,edenl,iondenl,etempl,eenrgl,ejl,elhl,inelhl))),delimiter="  ")
