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
        'planet0'       : pd.earth[ 'ID' ],     # Departure planet
        'planet1'       : pd.mars[ 'ID' ],      # Target planet
        'departure0'    : '2020-07-01',         # Intial departure date
        'departure1'    : '2020-09-01',         # Final departure date
        'arrival0'      : '2020-11-01',         # Initial arrival date
        'arrival1'      : '2022-01-24',         # Final arrival date
        'mu'            : pd.sun[ 'mu' ],       # Gravitational parameter in km**3/s**2
        'step'          : 1,                    # Step size in days
        'frame'         : 'J2000',              # Ecliptic of J2000
        'observer'      : '500@0',              # Solar Sytem Barycenter
        'cutoff_v'      : 20.0,                 # Maximum vinf to consider             
        'c3_levels'     : None,                 # C3 levels for contour plot
        'vinf_levels'   : None,                 # vinf levels for contour plot
        'tof_levels'    : None,                 # tof levels for contour plot
        'dv_levels'     : None,                 # dv levels for contour plot
        'dv_cmap'       : 'RdPu_r',             # color map for dv contours
        'figsize'       : ( 20, 10 ),           # figure size for contour plot
        'lw'            : 1.5,                  # linewidth for contour lines
        'title'         : 'Porkchop Plot',      # Plot title
        'fontsize'      : 15,                   # Axes fontsize
        'show'          : False,                # For displaying the figure
        'filename'      : None,                 # Specify filename for c3 plot
        'filename_dv'   : None,                 # Specify filename for dv plot
        'dpi'           : 300,                  # Specify target dpi
        'load'          : False                 #
    }
    
    # Overrides default config parameters
    for key in config.keys():
        _config[ key ] = config [ key ]

    # Define cutoff C3
    cutoff_c3 = _config[ 'cutoff_v' ] ** 2

    # Generate URLs for querying ephemeris data from Horizons API
    url_departure = generate_url(
        _config[ 'planet0' ],
        _config[ 'departure0' ],
        _config[ 'departure1' ],
        _config[ 'step' ]
    ) 
    url_arrival   = generate_url(
        _config[ 'planet1' ],
        _config[ 'arrival0' ],
        _config[ 'arrival1' ],
        _config[ 'step' ]
    )

    # Define target output path for departure and arrival data
    departure_output_path = f"{_config[ 'planet0' ]}_{_config[ 'departure0' ]}_{_config[ 'departure1']}.txt"
    arrival_output_path = f"{_config[ 'planet1' ]}_{_config[ 'arrival0' ]}_{_config[ 'arrival1']}.txt"

    # Submit API request and save the response to text files in target paths
    save_query_to_file(
        url_departure,
        departure_output_path
    )
    save_query_to_file(
        url_arrival,
        arrival_output_path
    )

    # Arrays of departure/arrival times and states
    et_departures, states_depart = stateReader(departure_output_path)
    et_arrivals, states_arrive   = stateReader(arrival_output_path)

    # Number of days in each array and total combinations
    ds  = len( et_departures )
    as_ = len( et_arrivals   )
    total = ds * as_

    print( 'Departure days: %i.'     % ds    )
    print( 'Arrival days: %i.'       % as_   )
    print( 'Total Combinations: %i.' % total )

    # Create empty array for C3, v_infinity, and tof to store data from lambert's solutions
    C3_shorts     = np.zeros( (as_, ds) )
    C3_longs      = np.zeros( (as_, ds) )
    v_inf_shorts  = np.zeros( (as_, ds) )
    v_inf_longs   = np.zeros( (as_, ds) )
    tofs          = np.zeros( (as_, ds) )

    # Create arrays for indexing and meshgrid
    x = np.arange( ds  )
    y = np.arange( as_ )

    # For each combination
    for na in y:
        for nd in x:

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