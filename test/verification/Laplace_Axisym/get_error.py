from sys import argv
import yt
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

ds=yt.load(argv[1])
prob_lo=ds.domain_left_edge.d
prob_hi=ds.domain_right_edge.d

data=np.loadtxt(argv[2])
x=data[:,0]
pot=data[:,1]

rmin=prob_lo[0]
rmax=prob_hi[0]
phi2=20.0
phi1=10.0

npts=len(x)
exactsoln=(phi2*np.log(x/rmin)+phi1*np.log(rmax/x))/np.log(rmax/rmin)
err_pot=np.sqrt(np.mean((pot-exactsoln)**2))
print("error:",err_pot)
plt.plot(x,exactsoln,label="exact",color="black")
plt.plot(x,pot,'r*',label="computed")
plt.legend()
plt.show()
