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

# Dark plotting background
plt.style.use( 'dark_background' )


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
        'step'          : 5,                    # Step size in days
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

    # Get ephemeris times and states
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

    # For each combination of departure and arrivals
    for nd, dep in enumerate( et_departures ):
        for na, arr in enumerate ( et_arrivals ):

            if arr > dep: # Ensure arrival date is after departure date
                # Calculate transfer time, seconds
                tof = (arr - dep) * 3600 * 24

            # Attempt to solve Lambert's problem for velocities

            # Short way (prograde)
            try:
                v_sc_depart_short, v_sc_arrive_short = lambert(
                    states_depart[ nd, :3 ],
                    states_arrive[ na, :3 ],
                    tof,
                    _config[ 'mu' ],
                    trajectory='pro'
                )
            except:
                v_sc_depart_short = np.array( [1000, 1000, 1000] )
                v_sc_arrive_short = np.array( [1000, 1000, 1000] )
            
            # Long way (retrograde)
            try:
                v_sc_depart_long, v_sc_arrive_long = lambert(
                    states_depart[ nd, :3 ],
                    states_arrive[ na, :3 ],
                    tof,
                    _config[ 'mu' ],
                    trajectory='retro'
                )
            except:
                v_sc_depart_long = np.array( [1000, 1000, 1000] )
                v_sc_arrive_long = np.array( [1000, 1000, 1000] )

            # Calculate C3 values at departure
            C3_short = norm( v_sc_depart_short - states_depart[ nd, 3: ] ) ** 2
            C3_long  = norm( v_sc_depart_long  - states_depart[ nd, 3: ] ) ** 2

            # Check for unreasonable values (C3)
            if C3_short > cutoff_c3: C3_short = cutoff_c3
            if C3_long  > cutoff_c3: C3_long = cutoff_c3

            # Calculate v_infinity values at arrival
            v_inf_short = norm( v_sc_arrive_short - states_arrive[ nd, 3: ] ) 
            v_inf_long  = norm( v_sc_arrive_long  - states_arrive[ nd, 3: ] )

            # Check for unreasonable values (v_infinity)
            if v_inf_short > _config[ 'cutoff_v' ]: v_inf_short = _config[ 'cutoff_v' ]
            if v_inf_long  > _config[ 'cutoff_v' ]: v_inf_long  = _config[ 'cutoff_v' ]

            # Append values to corresponding arrays
            C3_shorts    [ na, nd ] = C3_short
            C3_longs     [ na, nd ] = C3_long
            v_inf_shorts [ na, nd ] = v_inf_short
            v_inf_longs  [ na, nd ] = v_inf_long
            tofs         [ na, nd ] = tof

            print( f'{na} / {as_}.' )

    # Convert tof from sec to days
    tofs /= ( 3600.0 * 24.0 )

    # Total delta-v
    dv_shorts = v_inf_shorts + np.sqrt( C3_shorts )
    dv_longs  = v_inf_longs  + np.sqrt( C3_longs  )
    
    '''
    Plotting
    '''
    # Normalize the departure and arrival date grids 
    normed_departures = ( et_departures - et_departures[ 0 ] )
    normed_arrivals   = ( et_arrivals   - et_arrivals[ 0 ]   )

    # Generate departure and arrival date grids
    dep_mesh, arr_mesh = np.meshgrid( normed_departures, normed_arrivals )
    
    # Create levels arrays
    if _config  [ 'c3_levels'   ] is None:
        _config [ 'c3_levels'   ] = np.arange( 10, 50, 2)

    if _config  [ 'vinf_levels' ] is None:
        _config [ 'vinf_levels' ] = np.arange( 0, 15, 1)

    if _config  [ 'tof_levels'  ] is None:
        _config [ 'tof_levels'  ] = np.arange( 100, 500, 20)

    if _config  [ 'dv_levels'   ] is None:
        _config [ 'dv_levels'   ] = np.arange( 3, 20, 0.5)

    # Define linewdith
    lw = _config[ 'lw' ]

    # Create plots

    fig, ax = plt.subplots( figsize = _config[ 'figsize' ] )

    c0 = ax.contour(
        dep_mesh,
        arr_mesh,
        C3_shorts,
        levels = _config[ 'c3_levels' ], colors = 'm', linewidths = lw
    )
    c1 = ax.contour(
        dep_mesh,
        arr_mesh,
        C3_longs,
        levels = _config[ 'c3_levels' ], colors = 'm', linewidths = lw
    )
    c2 = ax.contour(
        dep_mesh,
        arr_mesh,
        v_inf_shorts,
        levels = _config[ 'vinf_levels' ], colors = 'deepskyblue', linewidths = lw
    )
    c3 = ax.contour(
        dep_mesh,
        arr_mesh,
        v_inf_longs,
        levels = _config[ 'vinf_levels' ], colors = 'deepskyblue', linewidths = lw
    )
    c4 = ax.contour(
        dep_mesh,
        arr_mesh,
        tofs,
        levels = _config[ 'tof_levels' ], colors = 'white', linewidths = lw * 0.6
    )

    
    plt.clabel( c0, fmt = '%i')
    plt.clabel( c1, fmt = '%i')
    plt.clabel( c2, fmt = '%i')
    plt.clabel( c3, fmt = '%i')
    plt.clabel( c4, fmt = '%i')

    plt.plot( [ 0 ], [ 0 ], 'm' )
    plt.plot( [ 0 ], [ 0 ], 'c' )
    plt.plot( [ 0 ], [ 0 ], 'w' )

    plt.legend(
        [
            r'C3 ($\dfrac{km^2}{s^2}$)',
            r'$V_{\infty}\; (\dfrac{km}{s})$',
            r'Time of Flight (days)'
        ],
        bbox_to_anchor = ( 1.005, 1.01 ),
        fontsize = 10
    )

    ax.set_title( _config[ 'title' ], fontsize = _config[ 'fontsize' ] )
    ax.set_ylabel( 'Arrival (Days Past %s)' % _config[ 'arrival0' ], fontsize = _config[ 'fontsize' ] )
    ax.set_xlabel( 'Departure (Days Past %s)' % _config[ 'departure0' ], fontsize = _config[ 'fontsize' ] )

    if _config[ 'show' ]:
        plt.show()

    if _config[ 'filename' ] is not None:
        plt.savefig( _config[ 'filename' ], dpi = _config[ 'dpi' ] )
        print( 'Saved', _config[ 'filename' ] )
    
    '''
    delta V plot
    '''

    pass