'''
Porkchop Plot Generator
'''

# Python Standard Libraries
import os

# 3rd Party Libraries
import numpy as np
import matplotlib.pyplot as plt

# Porkchop-Plot-Generator libraries
from utils import planetary_data  as pd
from utils import lambert_tools   as lt
from utils import ephemeris_query as eq
from utils.numerical_tools import norm

# Dark plotting background
plt.style.use( 'dark_background' )


def interplanetary_porkchop( departPlanet, targetPlanet, config ):
    
    # Default config dictionary
    _config = {
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
        'load'          : False                 # Load existing ephemeris data
    }
    
    # Overrides default config parameters
    for key in config.keys():
        _config[ key ] = config [ key ]

    '''
    Data handling and Ephemeris Query
    '''

    # Determine the directory for saving ephemeris data
    current_dir = os.path.dirname( __file__ )
    project_root = os.path.dirname( current_dir )
    data_dir = os.path.join( project_root, 'data' )

    # Create the main data directory if it doesn't exist
    if not os.path.exists( data_dir ):
        os.makedirs( data_dir, exist_ok = True )

    # Create subdirectories for departure and arrival data 
    departure_dir = os.path.join( data_dir, 'departure_data' )
    arrival_dir = os.path.join( data_dir, 'arrival_data')

    # Create the data subdirectories if they doesn't exist
    if not os.path.exists( departure_dir ):
        os.makedirs( departure_dir, exist_ok = True )

    if not os.path.exists( arrival_dir ):
        os.makedirs( arrival_dir, exist_ok = True )

    # Define target output path for departure and arrival data
    departure_output_path = os.path.join(
        departure_dir,
        f"{ departPlanet[ 'name'] }_{ _config[ 'departure0' ] }_{ _config[ 'departure1' ] }.txt"
    )
    arrival_output_path = os.path.join(
        arrival_dir,
        f"{ targetPlanet[ 'name' ] }_{ _config[ 'arrival0' ] }_{ _config[ 'arrival1' ] }.txt"
    )

    '''
    Check if the load parameter is set to True and if the necessary data files exist. If both conditions are met, it will load the data from these files instead of querying the API. If not, it will proceed with the API requests as usua
    '''
    if _config[ 'load' ] and os.path.exists( departure_output_path ) and os.path.exists( arrival_output_path ):
        print('Loading ephemeris data from existing files.')
    else:
        # Generate URLs for querying ephemeris data from Horizons API
        url_departure = eq.generate_url(
            departPlanet[ 'ID' ],
            _config[ 'departure0' ],
            _config[ 'departure1' ],
            _config[ 'step' ]
        ) 
        url_arrival   = eq.generate_url(
            targetPlanet[ 'ID' ],
            _config[ 'arrival0' ],
            _config[ 'arrival1' ],
            _config[ 'step' ]
        )

        # Submit API request and save the response to text files in target paths
        eq.save_query_to_file(
            url_departure,
            departure_output_path
        )
        eq.save_query_to_file(
            url_arrival,
            arrival_output_path
        )

    '''
    Calculations
    '''

    # Define cutoff C3
    cutoff_c3 = _config[ 'cutoff_v' ] ** 2

    # Get ephemeris times and states
    et_departures, states_depart = eq.stateReader(departure_output_path)
    et_arrivals, states_arrive   = eq.stateReader(arrival_output_path)

    # Number of days in each array and total combinations
    ds  = len( et_departures )
    as_ = len( et_arrivals   )
    total = ds * as_

    # Create empty array for C3, v_infinity, and tof to store data from lambert's solutions
    C3_shorts     = np.zeros( (as_, ds) )
    C3_longs      = np.zeros( (as_, ds) )
    v_inf_shorts  = np.zeros( (as_, ds) )
    v_inf_longs   = np.zeros( (as_, ds) )
    tofs          = np.zeros( (as_, ds) )

    # For each combination of departure and arrivals
    for na, arr in enumerate( et_arrivals ):
        for nd, dep in enumerate ( et_departures ):

            if arr > dep: # Ensure arrival date is after departure date
                # Calculate transfer time, seconds
                tof = (arr - dep) * 3600 * 24

            # Attempt to solve Lambert's problem for velocities

            # Short way (prograde)
            try:
                v_sc_depart_short, v_sc_arrive_short = lt.lambert_solver(
                    states_depart[ nd, :3 ],
                    states_arrive[ na, :3 ],
                    tof,
                    _config[ 'mu' ],
                    trajectory='pro'
                )
            except Exception as e:
                print(f"Prograde Lambert solution failed: {e}")
                v_sc_depart_short = np.array( [1000, 1000, 1000] )
                v_sc_arrive_short = np.array( [1000, 1000, 1000] )
            
            # Long way (retrograde)
            try:
                v_sc_depart_long, v_sc_arrive_long = lt.lambert_solver(
                    states_depart[ nd, :3 ],
                    states_arrive[ na, :3 ],
                    tof,
                    _config[ 'mu' ],
                    trajectory='retro'
                )
            except Exception as e:
                print(f"Retrograde Lambert solution failed: {e}")
                v_sc_depart_long = np.array( [1000, 1000, 1000] )
                v_sc_arrive_long = np.array( [1000, 1000, 1000] )

            # Calculate C3 values at departure
            C3_short = norm( v_sc_depart_short - states_depart[ nd, 3: ] ) ** 2
            C3_long  = norm( v_sc_depart_long  - states_depart[ nd, 3: ] ) ** 2

            # Check for unreasonable values (C3)
            if C3_short > cutoff_c3: C3_short = cutoff_c3
            if C3_long  > cutoff_c3: C3_long = cutoff_c3

            # Calculate v_infinity values at arrival
            v_inf_short = norm( v_sc_arrive_short - states_arrive[ na, 3: ] ) 
            v_inf_long  = norm( v_sc_arrive_long  - states_arrive[ na, 3: ] )

            # Check for unreasonable values (v_infinity)
            if v_inf_short > _config[ 'cutoff_v' ]: v_inf_short = _config[ 'cutoff_v' ]
            if v_inf_long  > _config[ 'cutoff_v' ]: v_inf_long  = _config[ 'cutoff_v' ]

            # Append values to corresponding arrays
            C3_shorts    [ na, nd ] = C3_short
            C3_longs     [ na, nd ] = C3_long
            v_inf_shorts [ na, nd ] = v_inf_short
            v_inf_longs  [ na, nd ] = v_inf_long
            tofs         [ na, nd ] = tof

        print( f'{(na + 1) / as_ * 100:.1f}%' )

    # After the loop, compute statistics for C3 and v_inf
    min_C3_short = np.min(C3_shorts)
    min_v_inf_short = np.min(v_inf_shorts)

    # Print results summary
    print("\n")
    print(f"{"*" * 51}")
    print(f"{" " * 18}Results Summary")          
    print(f"{"*" * 51}")
    print(f"Departure body name     : {departPlanet[ 'name' ]} ({departPlanet[ 'ID' ]})")
    print(f"Target body name        : {targetPlanet[ 'name' ]} ({targetPlanet[ 'ID' ]})")
    print(f"Center body name        : SOLAR SYSTEM BARYCENTER")
    print(f"Reference frame         : Ecliptic of {_config[ 'frame' ]}")
    print(f"{"*" * 51}")
    print(f"Launch window           : {_config[ 'departure0' ]} --> {_config[ 'departure1' ]}")
    print(f"Arrival window          : {_config[ 'arrival0' ]} --> {_config[ 'arrival1' ]}" )
    print(f"Step-size               : {_config[ 'step' ]}")
    print(f"{"*" * 51}")
    print(f'Departure days          : {ds}'    )
    print(f'Arrival days            : {as_}'   )
    print(f'Total Combinations      : {total}' )
    print(f'Trajectory Calculations : {2*total}')
    print(f"{"*" * 51}")
    print(f'Minimum departure C3    : {min_C3_short:.2f} (km**2/s**2)')
    print(f'Minimum arrival delta-V : {min_v_inf_short:.2f}  (km/s)')

    # Convert tof from sec to days
    tofs /= ( 3600.0 * 24.0 )

    # Total delta-v
    dv_shorts = v_inf_shorts + np.sqrt( C3_shorts )
    dv_longs  = v_inf_longs  + np.sqrt( C3_longs  )
    
    '''
    Plotting
    '''

    # Create subdirectories for figures
    fig_dir = os.path.join( data_dir, 'fig' )

    # Create the fig subdirectory if it doesn't exist
    if not os.path.exists( fig_dir ):
        os.makedirs( fig_dir, exist_ok = True )

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

    ax.set_title(
        f"{departPlanet[ 'name' ]} to {targetPlanet[ 'name' ]} Porkchop Plot",
        fontsize = _config[ 'fontsize' ]
    )
    ax.set_ylabel( 'Arrival (Days Past %s)' % _config[ 'arrival0' ], fontsize = _config[ 'fontsize' ] )
    ax.set_xlabel( 'Departure (Days Past %s)' % _config[ 'departure0' ], fontsize = _config[ 'fontsize' ] )

    if _config[ 'filename' ] is not None:
        plt.savefig( os.path.join( fig_dir, _config[ 'filename' ] ), dpi = _config[ 'dpi' ] )
        print( 'Saved', _config[ 'filename' ] )

    if _config[ 'show' ]:
        plt.show()
    
    plt.close()

    '''
    Total delta V plot
    '''

    fig, ax = plt.subplots( figsize = _config[ 'figsize' ] )

    c0 = ax.contour(
        dep_mesh,
        arr_mesh,
        dv_shorts,
        levels = _config[ 'dv_levels' ],
        cmap = _config[ 'dv_cmap' ], 
        linewidths = lw
    )
    c1 = ax.contour(
        dep_mesh,
        arr_mesh,
        dv_longs,
        levels = _config[ 'dv_levels' ],
        cmap = _config[ 'dv_cmap' ], 
        linewidths = lw
    )
    c2 = ax.contour(
        dep_mesh,
        arr_mesh,
        tofs,
        levels = _config[ 'tof_levels' ],
        colors = 'c', 
        linewidths = lw * 0.6
    )

    plt.clabel( c0, fmt = '%.1f' )
    plt.clabel( c1, fmt = '%.1f' )
    plt.clabel( c2, fmt = '%i' )

    ax.set_title(
        rf"{departPlanet[ 'name' ]} to {targetPlanet[ 'name' ]}: Total $\Delta V \; \left(\dfrac{{km}}{{s}}\right)$ Plot",
        fontsize = _config[ 'fontsize' ]
        )
    ax.set_ylabel( 'Arrival (Days Past %s)' % _config[ 'arrival0' ], fontsize = _config[ 'fontsize' ] )
    ax.set_xlabel( 'Departure (Days Past %s)' % _config[ 'departure0' ], fontsize = _config[ 'fontsize' ] )

    if _config[ 'filename_dv' ] is not None:
        plt.savefig( os.path.join( fig_dir, _config[ 'filename_dv' ] ), dpi = _config[ 'dpi' ] )
        print( 'Saved', _config[ 'filename_dv' ] )

    if _config[ 'show' ]:
        plt.show()

    plt.close()