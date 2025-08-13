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

outfile=open("avgdata","w")

for i, fn in enumerate(fn_list):

    ds=yt.load(fn)
    ad=ds.all_data()
    ref=np.array(ad["ReducedEF"])
    etemp=np.array(ad["Electron_Temp"])
    eden=np.array(ad["E"])
    Arpden=np.array(ad["ARp"])
    Ar2pden=np.array(ad["AR2p"])
    Hpden=np.array(ad["Hp"])
    H2pden=np.array(ad["H2p"])
    H3pden=np.array(ad["H3p"])
    ARHpden=np.array(ad["ARHp"])

    Hden=np.array(ad["H"])
    H2v1den=np.array(ad["H2v1"])
    H2v2den=np.array(ad["H2v2"])
    H2v3den=np.array(ad["H2v3"])
    
    outfile.write("%e\t%e\t%e\t%e\t"%(float(ds.current_time),np.mean(ref),np.mean(etemp),np.mean(eden)))
    outfile.write("%e\t%e\t%e\t%e\t"%(np.mean(Arpden),np.mean(Ar2pden),np.mean(Hpden),np.mean(H2pden)))
    outfile.write("%e\t%e\t%e\t%e\t"%(np.mean(H3pden),np.mean(ARHpden),np.mean(Hden),np.mean(H2v1den)))
    outfile.write("%e\t%e\n"%(np.mean(H2v2den),np.mean(H2v3den)))

outfile.close()
