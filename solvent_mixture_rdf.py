# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a function script file.
This file is used for RDF calculation

"""

import MDAnalysis as mda
import pandas as pd
import numpy as np
import os.path

import matplotlib.pyplot as plt
from gridData import Grid

# import nglview as nv
# %matplotlib inline

from MDAnalysis.analysis import rdf
import matplotlib.pyplot as plt

def system_sep():

    sys_co = mda.Universe(os.path.join("data_nvt_samp.lammpsdata"),
    os.path.join("dump_nvt_samp.lammpsdata.lammpstrj"), 
    topology_format="data", atom_style="id mol type q x y z", format="lammpsdump")

    zeo_atom_type = " ".join (str(i) for i in range (7,13))
    ad_atom_type= " ".join (str(i) for i in (range (13,24)))
    
    zeo = sys_co.select_atoms(f'type {ad_atom_type}')
    ad = sys_co.select_atoms(f'type {zeo_atom_type}')

    water = sys_co.select_atoms("type 2 5")
    methanol = sys_co.select_atoms("type 1 3 4 6")
    solvent = sys_co.select_atoms("type 1 2 3 4 5 6")
    
    return(zeo, ad, water, methanol, solvent)

def atom_rdf(atom1, atom2, start=None, stop=None, setstep=10, setnbin=100,range_ini=0.01):
    atom1_atom2_rdf = rdf.InterRDF(atom1, atom2,
                    nbins=setnbin,  # default
                    range=(range_ini, 15.0),  # distance
                    norm='density',
                              )
    
    atom1_atom2_rdf.run(start=None, stop=None,step=setstep)

#     plt.plot(atom1_a?tom2_rdf.bins, atom1_atom2_rdf.rdf)

    return (atom1_atom2_rdf.bins, atom1_atom2_rdf.rdf)

    
def atom_rdf_sol():
    
    system = system_sep()
    water = system[2]
    methanol = system[3]
    solvent = system[4]
    zeo = system[0]
    Owm=solvent.select_atoms("type 5 6" )
    Ow=water.select_atoms("type 5" )
    Om=methanol.select_atoms("type 6")
#     Ow=water.select_atoms(f'type 5 and sphzone {r_ti} type 7', updating=True)
#     Om=methanol.select_atoms(f'type 6 and sphzone {r_ti} type 7', updating=True)    
    
    Ti=zeo.select_atoms("type 7")
    
    
    
    return (Ow, Om, Ti, Owm)
#     return(system)


############################################ for pure water #################################33
def system_sep_water():
    sys_co = mda.Universe(os.path.join("data_nvt_samp.lammpsdata"),
    os.path.join("dump_nvt_samp.lammpsdata.lammpstrj"), 
    topology_format="data", atom_style="id mol type q x y z", format="lammpsdump")

    zeo_atom_type = " ".join (str(i) for i in range (3,9))
    ad_atom_type= " ".join (str(i) for i in (range (9,20)))
          
    zeo = sys_co.select_atoms(f'type {ad_atom_type}')
    ad = sys_co.select_atoms(f'type {zeo_atom_type}')
    water = sys_co.select_atoms("type 1 2")
    
    return(zeo, ad, water)

def atom_rdf_sol_water(r_ti=10):
    
    system = system_sep_water()
    water = system[2]
    zeo = system[0]
    Ow=water.select_atoms(f'type 1')
#     Ow=water.select_atoms(f'type 5 and sphzone {r_ti} type 7', updating=True)
#     Om=methanol.select_atoms(f'type 6 and sphzone {r_ti} type 7', updating=True)    
    
    Ti=zeo.select_atoms(f'type 3')
    
    return (Ow, Ti)


####################################

if __name__ == "__main__":
    ini_wd= r"C:\Users\chen.13058\OneDrive - The Ohio State University\Research_osu\project\zeolite\FAU_cosolvent"
    os.chdir(ini_wd)
    path_to_data=r"03_methanol_water_240_960\02_lammps_generate\10_free\01_methanol\411"
    os.chdir(os.path.join(ini_wd,path_to_data))
    system_sep()