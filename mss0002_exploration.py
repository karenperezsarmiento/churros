import pixell 
from pixell import enmap
from pixell import enplot
import numpy as np
import re

sim_dir = "/global/cfs/cdirs/sobs/sims/mss-0002/RC1.r01"
sim_fn = "sobs_RC1.r01_LAT_mission_%s_4way_%s_%s_map_%s.fits" #{(freq,splitnum,type,projection)}

out_dir = "/global/homes/k/kaper/mss02_exploration/"

projection = "car"
frequencies = ["f030","f040","f090","f150","f230","f290"]
types = ["signal","noise","sky_ivar","sky_map"]
splitnum = ["coadd","split1","split2","split3"]

sky_map_fn = f"{sim_dir}{sim_fn%(projection,frequencies[2],types[3],splitnum[0])}"
print(sky_map_fn)
sky_map = enmap.read_map(sky_map_fn)
print(sky_map.shape)
print(sky_map.wcs)
sky_d8 = enmap.downgrade(sky_map,8)
out_fn = re.sub(".fits","",sky_map_fn)
out_fn = out_fn + "_d8"
enplot.write(f"{out_dir}{out_fn}",enplot.plot(sky_d8,grid=True))



ivar_map_fn = f"{sim_dir}{sim_fn%(projection,frequencies[2],types[2],splitnum[0])}"
print(ivar_map_fn)
ivar_map = enmap.read_map(ivar_map_fn)
print(ivar_map.shape)
print(ivar_map.wcs)
ivar_d8 = enmap.downgrade(ivar_map,8,op=np.sum)
out_fn = re.sub(".fits","",ivar_map_fn)
out_fn = out_fn + "_d8"
enplot.write(f"{out_dir}{out_fn}",enplot.plot(sky_d8,grid=True))