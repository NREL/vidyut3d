import yt
from sys import argv
import matplotlib.pyplot as plt
import numpy as np
import glob

fn_pattern= argv[1]
fieldname=argv[2]
fn_list=[]
try:
    fn_list = sorted(glob.glob(fn_pattern), key=lambda f: int(f.split("plt")[1]))
except:
    if(fn_list==[]):
        print("using file of plotfiles..")
        infile=open(argv[1],'r')
        for line in infile:
            fn_list.append(line.split()[0])
        infile.close()

first_res=1.0
resarr=np.zeros((len(fn_list)-1,4))
for i in range(1,len(fn_list)):

    dsm1=yt.load(fn_list[i-1])
    ds=yt.load(fn_list[i])
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
    fields_load=[fieldname]
    ad   = ds.covering_grid(level=covgrid_lev, left_edge=prob_lo, dims=res, fields=fields_load)
    adm1 = dsm1.covering_grid(level=covgrid_lev, left_edge=prob_lo, dims=res, fields=fields_load)

    field=np.array(ad[fieldname])
    fieldm1=np.array(adm1[fieldname])
    residual=np.abs(np.max(field)-np.max(fieldm1))
    if(i==1):
        first_res=residual
    print("abs/rel residual:",residual,residual/first_res)
    resarr[i-1,0]=i
    resarr[i-1,1]=residual
    resarr[i-1,2]=residual/first_res
    resarr[i-1,3]=np.mean(field)


np.savetxt("resnorms",resarr,delimiter="  ")
