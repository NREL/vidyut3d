import yt
from sys import argv
import matplotlib.pyplot as plt
import numpy as np
import glob
    
def _efieldmag(field, data):
    return (np.sqrt(
        data["Efieldx"]*data["Efieldx"]+ 
        data["Efieldy"]*data["Efieldy"]+
        data["Efieldz"]*data["Efieldz"])); 
    
fn_pattern= argv[1]
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

fieldname="efieldmag"
set_minmax=False

print(fn_list)
if(len(argv) > 3):
    minval=float(argv[2])
    maxval=float(argv[3])
    set_minmax=True

for i, fn in enumerate(fn_list):
    ds=yt.load(fn)
    ds.add_field(("gas", "efieldmag"), function=_efieldmag, units="", sampling_type="local")
    slc = yt.SlicePlot(ds, 'z', fieldname)
    slc.set_log(fieldname,False)
    if(set_minmax):
        slc.set_zlim(fieldname,minval,maxval)
    slc.annotate_grids()
    slc.save(fieldname+"_%4.4d"%(i))
