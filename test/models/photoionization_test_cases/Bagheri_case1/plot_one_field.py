import numpy as np
from matplotlib.colors import LogNorm
import os 
import fnmatch
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker,cm

import yt
from yt.visualization.base_plot_types import get_multi_plot


orient = "vertical"

plot_size = 0.0125
resolution = 1.5e-6 # [m]

nx = int(plot_size/resolution)
ny = nx

# Select the time 
istart = 0 # Start index
iskip  = 1 # Skip index
path = os.getcwd()
# path = '/lustre/scratch/ldowen/Bruno/'

list_output=fnmatch.filter(os.listdir(path), 'plt*')
list_output.sort()    
print("Found output   : {}".format(list_output))

list_process=list_output[istart:len(list_output):iskip]
list_process=[list_process[-1]]
#list_process=list_process[-2:]
list_process=['plt00011']
print("To be processed: {} \n".format(list_process))


scalars = ["E"]
#scalars = ['Electron_Temp','N2p','ReducedEF']

for i,name in enumerate(list_process):
        # ds = yt.load("{}/{}".format(path,name))


    ds1 = yt.load("{}/{}".format(path,name))  # load data
    time = float(ds1.current_time)
    # scalars = [name[1] for name in ds1.field_list]
    ds1.print_stats()

    for scalar in scalars:

        slc1 = yt.SlicePlot(
            ds1,
            "theta",
            #center=[1.0e-4,0.0,30.e-3],
            center = "c",
            fields=[scalar],
        )

        #slc1_Z = yt.SlicePlot(
        #    ds1,
        #    "z",
        #    center=[0.0,0.0,0.0],
        #    #center = "c",
        #    fields=[scalar],
        #)

        slc_frb1 = slc1.data_source.to_frb((plot_size), nx)
        #slc_frb1_z = slc1_Z.data_source.to_frb((plot_size), nx)


        slc_temp1 = np.array(slc_frb1[scalar])
        #slc_tempz = np.array(slc_frb1_z[scalar])
        #slc_vfrac = np.array(slc_frb1["volFrac"])
        
        #I = np.where(slc_frb1["volFrac"] < 0.01)
        #slc_temp1[I] = np.nan 
     
        extent = -plot_size*1.e+3/2.0, plot_size*1.e+3/2.0, -10.0, (plot_size)*1.e+3-10.0
        extentZ = -plot_size*1.e+3/2.0, plot_size*1.e+3/2.0, -plot_size*1.e+3/2.0, plot_size*1.e+3/2.0
        xgrid = np.linspace( -plot_size*1.e+3/2.0, plot_size*1.e+3/2.0, nx)
        ygrid = np.linspace( 0.0                 , plot_size*1.e+3    , ny)

        fig, ax1 = plt.subplots()
        fig2, ax2 = plt.subplots()
        minval = np.nanmin(slc_temp1)
        maxval = np.nanmax(slc_temp1)
        # CS = ax1.contourf(xgrid,ygrid,slc_temp1,levels=np.linspace(0.,8.e-3,20),cmap=cm.turbo)
        CS = ax1.imshow(slc_temp1,vmin=minval,vmax=maxval,cmap=cm.turbo,extent=extent,origin='lower')
        #CSz = ax2.imshow(slc_tempz,vmin=minval,vmax=maxval,cmap=cm.turbo,extent=extentZ,origin='lower')
        cbar = fig.colorbar(CS, ax=ax1, format='%03.1e', ticks=np.linspace(minval,maxval,10),shrink=0.8)
        cbar.ax.set_ylabel(scalar)

        cbar = fig.colorbar(CS, ax=ax2, format='%03.1e', ticks=np.linspace(minval,maxval,10),shrink=0.8)
        cbar.ax.set_ylabel(scalar)

        ax1.set_xlabel('x [mm]')
        ax1.set_ylabel('z [mm]')
        ax1.set_xlim(-plot_size*1.e+3/2.0,plot_size*1.e+3/2.0)
        ax1.set_ylim(-10.,plot_size*1.e+3-10.)

        fig.savefig('{}_{:1.2f}ns.png'.format(scalar,time*1.0e+9),dpi=300)

plt.show()
