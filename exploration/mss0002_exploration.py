import pixell 
from pixell import enmap
from pixell import enplot
import numpy as np
import re

sim_dir = "/gpfs/fs0/project/r/rbond/kaper/so/mss0002/mss0002_sims/"
sim_fn = "sobs_RC1.r01_LAT_mission_%s_4way_%s_%s_%s.fits" #{(freq,splitnum,type,projection)}

out_dir = "/home/s/sievers/kaper/scratch/mss002/plots/"

projection = "car"
frequencies = ["f030","f040","f090","f150","f230","f290"]
types = ["signal_map","noise_map","sky_ivar","sky_map"]
splitnum = ["coadd","split1","split2","split3"]

sky_map_fn = f"{sim_fn%(frequencies[2],splitnum[0],types[3],projection)}"
print(sky_map_fn)
sky_map = enmap.read_map(f"{sim_dir}{sky_map_fn}")
print(sky_map.shape)
print(sky_map.wcs)
sky_d8 = enmap.downgrade(sky_map,8)
out_fn = re.sub(".fits","",sky_map_fn)
out_fn = out_fn + "_d8"
enplot.write(f"{out_dir}{out_fn}",enplot.plot(sky_d8,grid=True,ticks=10))



ivar_map_fn = f"{sim_fn%(frequencies[2],splitnum[0],types[2],projection)}"
print(ivar_map_fn)
ivar_map = enmap.read_map(f"{sim_dir}{ivar_map_fn}")
print(ivar_map.shape)
print(ivar_map.wcs)
ivar_d8 = enmap.downgrade(ivar_map,8,op=np.sum)
out_fn = re.sub(".fits","",ivar_map_fn)
out_fn = out_fn + "_d8"
enplot.write(f"{out_dir}{out_fn}",enplot.plot(sky_d8,grid=True,ticks=10))