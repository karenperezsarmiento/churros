"""
Enhance masks with source cuts
"""


# %%
from pixell import enmap, curvedsky as cs, utils as u,reproject,wcsutils
from orphics import maps,io,cosmology,stats,pixcov
import numpy as np
import os,sys,warnings
import healpy as hp

# ARGPARSE

import argparse
# Parse command line
parser = argparse.ArgumentParser(description='Make a mask.')
parser.add_argument("outname", type=str,help='Name of outputs to enhance. Could include a path.')
parser.add_argument("--width-deg",     type=float,  default=3.0,help="Width in deg. to grow mask by.")
parser.add_argument("--planck-cuts",     type=str,  default="60,70,80",help="Galactic cuts, comma spearated.")
parser.add_argument("--decmin",     type=float,  default=None,help="Zero out all regions with dec less than this")
parser.add_argument("--decmax",     type=float,  default=None,help="Zero out all regions with dec more than this")
args = parser.parse_args()

cuts = args.planck_cuts.split(',')
fname = f"{args.outname}_joint_mask_{cuts[0]}.fits"
shape,wcs = enmap.read_map_geometry(fname)
shape = shape[-2:]
mask = enmap.ones(shape,wcs)

m1 = maps.mask_srcs(shape,wcs,np.asarray([[-5],[79]]),6.*60)
m2 = maps.mask_srcs(shape,wcs,np.asarray([[-5],[-121]]),3.5*60)
m3 = maps.mask_srcs(shape,wcs,np.asarray([[-13.0667],[-61.733]]),4.*60)
m4 = maps.mask_srcs(shape,wcs,np.asarray([[17.15],[108.633]]),5.*60)
m5 = maps.mask_srcs(shape,wcs,np.asarray([[1.35],[115.2]]),3.*60)
            
jmask = mask.astype(bool) * m1.astype(bool) * m2.astype(bool) * m3.astype(bool) * m4.astype(bool) * m5.astype(bool)
print("Applying large circle masks...")
gmask = maps.grow_mask(jmask,5)
omask = 1-maps.grow_mask(1-gmask.astype(int),5)

io.hplot(omask,f"{args.outname}_mask_enhancement1",downgrade=16,colorbar=True,ticks=2,grid=True,font_size=6)
p545 = enmap.ones(shape,wcs)

print("Applying 545 GHz cut in problem area...")
coords = np.asarray([
[-74.5,-36.9],
])[:,[1,0]] * u.degree

radius = 10 * u.degree

gals = pixcov.extract_cutouts(p545,coords,radius)


carray = np.asarray([
[-76,-38,0.5],
[-75,-37.2,1.1],
[-72,-37.4,1.5],
[-70,-37.2,2.4],
[-68,-36.2,1.8],
[-65.5,-36.2,2.],
[-63.0,-36.5,0.5],
[-64.0,-34.5,0.5],


])
def mask_circles(shape,wcs,carray):
    mask = enmap.ones(shape,wcs,dtype=bool)
    for c in carray:
        mask = mask * maps.mask_srcs(mask.shape,mask.wcs,c[:2][::-1][...,None],c[2]*60)
    return mask

for gal in gals:
    g = gal.copy()
    m = mask_circles(g.shape,g.wcs,carray)
    m = maps.grow_mask(m,2)
    m = 1 - maps.grow_mask(1-m,2)
    m = m.astype(float)
    

gmask = omask.copy()*0 + 1.
fmask = enmap.insert(gmask,m) * omask

def fsky(mask):
    return (mask.pixsize()*mask).sum() / np.pi / 4.

print("Enhancing and writing masks...")

for cut in cuts:
    gmap = enmap.read_map(f"{args.outname}_joint_mask_{cut}.fits")
    gmap[fmask<=0.5] = 0
    if args.decmin or args.decmax:
        dec_map, ra_map = enmap.posmap(gmap.shape,gmap.wcs)
    if args.decmin:
        gmap[dec_map < (args.decmin*u.degree)] = 0
    if args.decmax:
        gmap[dec_map > (args.decmax*u.degree)] = 0

    gmap = maps.grow_mask(gmap,args.width_deg)
    gmap = 1-maps.grow_mask(1-gmap.astype(int),args.width_deg)
    
    print(cut, fsky(gmap))
    enmap.write_map(f"{args.outname}_enhanced_mask_{cut}.fits",gmap)
    io.hplot(gmap,f"{args.outname}_enhanced_mask_{cut}",downgrade=16,colorbar=True,ticks=2,grid=True,font_size=6)