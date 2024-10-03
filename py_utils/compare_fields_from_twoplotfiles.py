import yt
import sys
from sys import argv
import matplotlib.pyplot as plt
import numpy as np
import glob

try:
    fn1= argv[1]
    fn2= argv[2]
    fieldname1=argv[3]
    fieldname2=argv[4]
except:
    print("\n\n*****************************************************\n")
    print("This script compares field variables from two plot files \nand prints the error between the data values")
    print("way to run:") 
    print("python compare_field_from_twoplotfiles.py <plotfile1> <plotfile2> <varname1> <varname2>")
    print("e.g. python compare_field_from_twoplotfiles.py plt1 plt2 Potential Potential")
    print("\n*****************************************************\n")
    sys.exit()

ds1=yt.load(fn1)
ds2=yt.load(fn2)
#both files should have same dimensions
prob_lo=ds1.domain_left_edge.d
prob_hi=ds1.domain_right_edge.d
probsize=prob_hi-prob_lo

ncells=ds1.domain_dimensions
    
minlev=0
maxlev=ds1.index.max_level
lengths=prob_hi-prob_lo
covgrid_lev=maxlev
res=np.array([ncells[0]* (2**covgrid_lev),ncells[1]* (2**covgrid_lev),ncells[2]* (2**covgrid_lev)])
dx_frb=probsize/res
fields_load1=[fieldname1]
fields_load2=[fieldname2]
ad1 = ds1.covering_grid(level=covgrid_lev, left_edge=prob_lo, dims=res, fields=fields_load1)
ad2 = ds2.covering_grid(level=covgrid_lev, left_edge=prob_lo, dims=res, fields=fields_load2)

field1=np.array(ad1[fieldname1])
field2=np.array(ad2[fieldname2])
err=np.sqrt(np.mean((field1-field2)**2))
field1avg=np.mean(np.abs(field1))
print("absolute and relative error for fields ("+fieldname1+"  "+fieldname2+") = %e\t%e"%(err,err/field1avg))
