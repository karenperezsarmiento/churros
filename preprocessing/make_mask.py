"""
                                                                   %@@.      :@@.                   
                                                                #@.              #@                 
                                                              -@                   #%               
                                                             @=     :@         .%:   @              
                                                            @      =@@:        =@@:   @             
                                                          :@              :-          ++            
                                                         @.              :@@+          @            
                                                      +@.             % :  @  = :      @            
                                                   =@+                =@@@@%@%%@       @            
                                               :@@.                                   *=    @       
                                           =@@.                                       @   =@ @      
                                       @@#                                           @  %@  @       
                                   *@*               @                              **:    @.       
                                @@.                  +=                            =%    :@         
                            .@%                       @.                          =%   #@           
                          @@                    @@%    %#                        @@@@=              
                       @@                          :@@@@                        @                   
                    :@:                                                       *@                    
                  .@                                                        :@                      
     @#.   .=@@= #*                                                       =@.                       
     :@@@%      :*                                                      @%                          
      @@%                                                           .@@                             
    -#                                                         .@@@                                 
     @@   :=#%%%%%%*+-.                          :-*%%@@@@@#:                                       

Lensing analysis mask preparation

To make a mask, we need the following things:
1. A Planck Galaxy mask
2. The ivar maps of the ACT/SO arrays
3. Lists of annoying regions to mask. This will be handled in a separate script.

This does step 2
"""


from pixell import enmap, curvedsky as cs, utils as u,reproject,wcsutils
import numpy as np
import os,sys,warnings,glob
import healpy as hp
from orphics import maps, io
from concurrent import futures


# ARGPARSE

import argparse
# Parse command line
parser = argparse.ArgumentParser(description='Make a mask.')
parser.add_argument("outname", type=str,help='Name of outputs. Could include a path.')
parser.add_argument("--rms-threshold",     type=float,  default=70.0,help="RMS threshold in uK-arcmin for ivar maps.")
parser.add_argument("--order",     type=int,  default=0,help="Planck mask interpolation order.")
parser.add_argument("--width-deg",     type=float,  default=0.2,help="Width in deg. to grow mask by.")
parser.add_argument("--nworkers",     type=int,  default=None,help="Maximum number of workers to parallelize over. Defaults to number of jobs.")
parser.add_argument("--template-fname",     type=str,  default=None,help="Path to template ACT/SO map.",required=True)
parser.add_argument("--planck-base-name",     type=str,  default=None,help="Root file name to Planck cut maps",required=True)
parser.add_argument("--planck-cuts",     type=str,  default="60,70,80",help="Galactic cuts, comma spearated.")
parser.add_argument("--ivar-search-string",     type=str,  default="cmb_night_?_3pass_4way_*_ivar.fits",help="Search string to find inverse variance maps in args.exp-path. The question mark will be replaced by the items in args.arrays.")
parser.add_argument("--arrays",     type=str,  default='pa4_f220,pa4_f150,pa5_f090,pa5_f150,pa6_f090,pa6_f150',help="Comma separated list of array names to replace * in args.ivar_search_string.")
parser.add_argument("--exp-path",     type=str, help="Path to ACT/SO maps.",required=True)
args = parser.parse_args()

jobs = []

# Get shape,wcs of output mask(s)
shape,wcs = enmap.read_map_geometry(args.template_fname)
owcs = wcs
shape = shape[-2:]

def load_map(i):
    desc = jobs[i]
    print("Starting ", desc)
    fname = desc[1]
    ivar = enmap.read_map(fname)
    if ivar.ndim==3:
        ivar = ivar[0]
    if ivar.ndim!=2: raise ValueError
    if ivar.shape[0]!=shape[0] or ivar.shape[1]!=shape[1]: raise ValueError(f"{fname} has shape {ivar.shape} but template has shape {shape}")
    if not(wcsutils.equal(ivar.wcs,owcs)): raise ValueError
    out = maps.rms_from_ivar(ivar)
    out[ivar<=0] = np.inf
        
    print("Done with ", desc)
    return out

cuts = args.planck_cuts.split(',')
if len(cuts)==1 and cuts[0]=='': cuts = []

for cut in cuts:
    if not(os.path.exists(args.planck_base_name+f"_{cut}.fits")): raise FileNotFoundError


arrays = args.arrays.split(',')
ivar_searches = [os.path.join(args.exp_path, args.ivar_search_string.replace('?',a)) for a in arrays]
ivar_files = []
for ivs in ivar_searches:
    ivar_files = ivar_files + glob.glob(ivs)
    
if len(ivar_files)<1:
    warnings.warn("No ivar files found.")
else:
    print("Will calculate RMS from the following ivar_files: ", ivar_files)

for ivar_file in ivar_files:
    if not(os.path.exists(ivar_file)): raise FileNotFoundError
    jobs.append(['rms',ivar_file])
    
njobs = len(jobs)

if args.nworkers is None:
    warnings.warn(f"Using as many workers as there are jobs: {njobs}. You may run out of memory.")

with futures.ProcessPoolExecutor(max_workers=args.nworkers) as executor:
    omaps = list(executor.map(load_map, range(njobs)))
    executor.shutdown(wait=True)



imask = enmap.ones(shape,wcs,dtype=int)

for i,desc in enumerate(jobs):
    if desc[0]=='rms':
        rms = omaps[i]
        imask[rms>=args.rms_threshold] = 0
enmap.write_map(f"{args.outname}_rms_mask.fits",imask)
io.hplot(imask,f"{args.outname}_rms_mask",downgrade=16,colorbar=True)


for cut in cuts:
    gal = enmap.extract(enmap.read_map(args.planck_base_name+f"_{cut}.fits"),imask.shape,imask.wcs)
    gal[imask==0] = 0.
    smask = maps.grow_mask(gal,args.width_deg)
    enmap.write_map(f"{args.outname}_joint_mask_{cut}.fits",smask.astype(int))
    io.hplot(smask,f"{args.outname}_joint_mask_{cut}",downgrade=16,colorbar=True)

        
        
    
    





