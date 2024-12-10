# Taaresh - Plots contours of axisymmetric data and mid-planes of 3D data

import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker,cm
import yt
from yt.visualization.base_plot_types import get_multi_plot
import os 
import fnmatch


orient = "vertical"

plot_size = 0.0125
resolution = 1.5e-6 # [m]

nx = int(plot_size/resolution)
ny = nx

# Select the time 
istart = 0 # Start index
iskip  = 1 # Skip index
path = os.getcwd()

list_output=fnmatch.filter(os.listdir(path), 'plt*')
list_output.sort()    
print("Found output   : {}".format(list_output))

list_process=list_output[istart:len(list_output):iskip]
list_process=[list_process[-1]]
#list_process=list_process[-2:]
list_process=['plt00075']
print("To be processed: {} \n".format(list_process))

#scalars = ["ReducedEF"]
scalars = ['Efieldx','Efieldy','PhotoIon_Src','E','Electron_Temp','N2p','ReducedEF']

for i,name in enumerate(list_process):

    ds1 = yt.load("{}/{}".format(path,name))  # load data
    time = float(ds1.current_time)
    # scalars = [name[1] for name in ds1.field_list]
    ds1.print_stats()

    for scalar in scalars:

        slc1 = yt.SlicePlot(
            ds1,
            "theta",
            center = "c",
            fields=[scalar],
        )

        slc_frb1 = slc1.data_source.to_frb((plot_size), nx)
        slc_temp1 = np.array(slc_frb1[scalar])
     
        extent = 0.0, plot_size*1.e+3, -10.0, (plot_size)*1.e+3-10.0
        xgrid = np.linspace( 0.0, plot_size*1.e+3, nx)
        ygrid = np.linspace( 0.0, plot_size*1.e+3, ny)

        fig, ax1 = plt.subplots()
        minval = np.nanmin(slc_temp1)
        maxval = np.nanmax(slc_temp1)
        CS = ax1.imshow(slc_temp1,vmin=minval,vmax=maxval,cmap=cm.turbo,extent=extent,origin='lower')
        # CSlog = ax1.imshow(slc_temp1,norm=LogNorm(vmin=max(1e20, 1e-10), vmax=1e26),cmap=cm.rainbow,extent=extent,origin='lower')
        
        cbar = fig.colorbar(CS, ax=ax1, format='%03.1e', ticks=np.linspace(minval,maxval,10),shrink=0.8)
        # cbarlog = fig.colorbar(CSlog, ax=ax1, format='%03.1e', ticks=np.logspace(np.log10(max(1e20, 1e-10)), np.log10(1e26), num=9), shrink=0.8)
        cbar.ax.set_ylabel(scalar)
        # cbarlog.ax.set_ylabel(scalar)

        ax1.set_xlabel('r [mm]')
        ax1.set_ylabel('z [mm]')
        ax1.set_xlim(0.0,plot_size*1.e+3)
        ax1.set_ylim(-10.,plot_size*1.e+3-10.)

        fig.savefig('{}_{:1.2f}ns.png'.format(scalar,time*1.0e+9),dpi=300)

plt.show()
