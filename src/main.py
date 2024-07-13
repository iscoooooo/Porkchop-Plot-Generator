# Porkchop-Plot-Generator Libraries
from utils import planetary_data as pd
from porkchop import interplanetary_porkchop

# Main script

def main():

    # config parameters for porkchop plot generator
    config = {
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
        'figsize'       : ( 6, 10 ),            # figure size for contour plot
        'lw'            : 1.5,                  # linewidth for contour lines
        'title'         : 'Porkchop Plot',      # Plot title
        'fontsize'      : 15,                   # Axes fontsize
        'show'          : True,                 # For displaying the figure
        'filename'      : None,                 # Specify filename for c3 plot
        'filename_dv'   : None,                 # Specify filename for dv plot
        'dpi'           : 300,                  # Specify target dpi
        'load'          : False                 # Load existing ephemeris data
    }

    # Call porkchop plot generator
    interplanetary_porkchop( config )

if __name__ == "__main__":
    main()