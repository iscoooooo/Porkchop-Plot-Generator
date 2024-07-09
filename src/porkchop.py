'''
Porkchop Plot Generator
'''

# Python Standard Libraries
import numpy as np
import matplotlib.pyplot as plt

# Porkchop-Plot-Generator libraries
import planetary_data as pd
from ephemeris_query import *
from utils import *


def interplanetary_porkchop( config ):
    
    # Default config dictionary
    _config = {
        'planet0'       : pd.earth[ 'ID' ],
        'planet1'       : pd.mars[ 'ID' ],
        'departure0'    : '2020-07-01',
        'departure1'    : '2020-09-01',
        'arrival0'      : '2020-11-01',
        'arrival1'      : '2022-01-24',
        'mu'            : pd.sun[ 'mu' ],
        'step'          : 1,
        'frame'         : 'J2000',
        'observer'      : '500@0',
        'cutoff_v'      : 20.0,
        'c3_levels'     : None,
        'vinf_levels'   : None,
        'tof_levels'    : None,
        'dv_levels'     : None,
        'dv_cmap'       : 'RdPu_r',
        'figsize'       : ( 20, 10 ),
        'lw'            : 1.5,
        'title'         : 'Porkchop Plot',
        'fontsize'      : 15,
        'show'          : False,
        'filename'      : None,
        'filename_dv'   : None,
        'dpi'           : 300,
        'load'          : False
    }
    
    # Overrides default config parameters
    for key in config.keys():
        _config[ key ] = config [ key ]

    # Compute cutoff C3
    cutoff_c3 = _config[ 'cutoff_v' ] ** 2

    # Arrays of departure and arrival times


    # Number of days in each array and total combinations


    # Create empty array for C3, v_infinity, and tof


    # Create arrays for indexing and meshgrid


    # For each combination



            # State of planet0 at departure


            # State of planet1 at arrival


            # Calculate transfer time


            # Attempt to solve Lambert's problem for velocities


            # Calculate C3 values at departure


            # Check for unreasonable  (C3)


            # Calculate v_infinity values at arrival


            # Check for unreasonable values (v_infinity)


            # Append values to corresponding array


    # Convert tof from sec to days


    # Total delta-v

    #------------------------------------------------------
    '''
    Plotting
    '''
    
    # Create levels arrays


    pass