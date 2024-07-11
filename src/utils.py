'''
Utility functions
'''

import csv
import numpy as np

def stateReader(filespec):
    '''
    Reads ephemeris data from a text file in CSV format and extracts Julian dates and state vectors.

    The function expects the file to contain ephemeris data in a specific format.
    It looks for the start ("SOE") and end ("EOE") markers to determine the data range
    and extracts the data between these markers. The extracted data is then converted
    to numpy arrays for further processing.

    Parameters:
    filespec : str
        The file path to the txt file containing the ephemeris data.

    Returns:
    julianDates : ndarray
        An array of Julian dates.
    states : ndarray
        An array of state vectors, where each state vector contains position and velocity
        components (x, y, z, vx, vy, vz).
    '''

    # Initialize lists to store Julian Dates and states
    julianDates = []
    states = []

    expr1 = 'SOE' # text to look for in the file
    expr2 = 'EOE' # text to look for in the file

    # Loop to find where the data range is at
    N = 0 # line count
    idx1 = None
    idx2 = None

    with open(filespec,'r') as file:
        while True:
            dataline = file.readline() # read one line at a time

            N = N + 1 # find the line number

            if isinstance(dataline, str): # if the line is a string

                # found the 1st string that I'm seeking
                if expr1 in dataline:
                    idx1 = N
                
                # found the 2nd string that I'm seeking
                if expr2 in dataline:
                    idx2 = N

            # break loop if both indices are found
            if (idx1 is not None) and (idx2 is not None):
                break
        
        # scan through data, skipping header, and collecting data
        with open(filespec, 'r') as file:
            reader = csv.reader(file,delimiter=',')

            # skip header lines to idx1
            for _ in range(idx1):
                next(reader)

            # collect data starting from line idx1 to idx2
            for row in reader:
                if reader.line_num >= idx2:
                    break

                julianDates.append(float(row[0]))
                states.append([float(value) for value in row[2:8]])
            
        # Convert lists to numpy arrays for further processing
        julianDates = np.array(julianDates)
        states = np.array(states)

    return julianDates, states


def lambert(R1, R2, dt, mu, tol=1e-6, maxiter=10000, trajectory='pro'):
    '''
    Solves Lambert's problem to find the initial and final velocity vectors (V1, V2)
    for a given initial and final position vectors (R1, R2) and time of flight (dt).
    
    This Lambert Solver is based on Algorithm 5.2 from "Orbital Mechanics for Engineering Students,"
    Fourth Edition by Howard D. Curtis.
    
    Parameters:
    R1, R2 : ndarray
        Initial and final position vectors (km)
    dt : float
        Time of flight from R1 to R2 (s)
    mu : float
        Gravitational parameter (km^3/s^2)
    tol : int, optional
        Tolerance for the solver
    maxiter : int, optional
        Maximum number of iterations for Newton's method
    string : str, optional
        'pro' for prograde orbit, 'retro' for retrograde orbit (default is 'pro')
    
    Returns:
    V1, V2 : ndarray
        Initial and final velocity vectors (km/s)
    '''

    # Ensure tolerance and maximum iterations are of the correct type
    if not isinstance(tol, (float, int)) or tol <= 0:
        raise ValueError("Tolerance 'tol' must be a positive number.")
    
    if not isinstance(maxiter, int) or maxiter <= 0:
        raise ValueError("Maximum iterations 'maxiter' must be a positive integer.")

    ## Helper functions:
    
    def y(z):
        return r1 + r2 + A * (z * S(z) - 1) / np.sqrt(C(z))

    def F(z):
        return (y(z) / C(z)) ** 1.5 * S(z) + A *  np.sqrt(y(z)) - np.sqrt(mu) * dt

    def dFdz(z):
        if z == 0:
            return np.sqrt(2) / 40 * y(0) ** 1.5 + A / 8 * (np.sqrt(y(0)) + A * np.sqrt(1 / 2 / y(0)))
        else:
            return (y(z) / C(z)) ** 1.5 * (1 / 2 / z * (C(z) - 3 * S(z) / 2 / C(z)) + 3 * S(z) ** 2 / 4 / C(z)) + A / 8 * (3 * S(z) / C(z) * np.sqrt(y(z)) + A * np.sqrt(C(z) / y(z)))
        

    # Magnitudes of R1 and R2
    r1 = np.linalg.norm(R1)
    r2 = np.linalg.norm(R2)

    # Compute the cross product between R1 and R2
    cross12 = np.cross(R1, R2)

    # Compute the change in true anomaly
    dtheta = np.arccos(np.dot(R1, R2) / (r1 * r2))

    if trajectory == 'pro':  # For prograde trajectory
        if cross12[2] < 0:
            dtheta = 2 * np.pi - dtheta
    elif trajectory == 'retro':  # For retrograde trajectory
        if cross12[2] >= 0:
            dtheta = 2 * np.pi - dtheta
    else:
        raise ValueError("Please indicate whether trajectory is prograde or retrograde with keyword argument.")

    # Compute the auxiliary function, A(r1,r2,dtheta)
    A = np.sin(dtheta) * np.sqrt(r1 * r2 / (1 - np.cos(dtheta)))

    # Determine approximately where F(z,dt) changes sign, and
    # use that value of z as the starting value for z:
    z = 0.1  # Start with a small positive value
    while F(z) < 0:
        z += 0.1
        if z > 1e6:  # Prevent infinite loop in case of an issue
            raise ValueError("Failed to find a suitable starting z value")

    # Find z by iterating using Newton's method until convergence within the error tolerance:
    z = newtonRaphson(F, dFdz, z, tol, maxiter)

    # Check the solution
    if np.isnan(z) or np.isinf(z):
        solved = False
    else:
        solved = True 
    
    if solved:
        # Compute the Lagrangian coefficients
        f = 1 - y(z) / r1
        fdot = (np.sqrt(mu) / (r1 * r2)) * np.sqrt(y(z) / C(z)) * (z * S(z) - 1)
        g = A * np.sqrt(y(z) / mu)
        gdot = 1 - y(z) / r2

        # Compute the velocities V1 & V2
        V1 = 1 / g * (R2 - f * R1)
        V2 = 1 / g * (gdot * R2 - R1)
            
        return V1, V2
    else: 
        print('Lambert solver did not converge.')
        return np.array([0, 0, 0]), np.array([0, 0, 0])

def C(z):
    '''
    Stumpff Function
    '''

    if z > 0:
        return (1 - np.cos(np.sqrt(z))) / z
    elif z < 0:
        return (np.cosh(np.sqrt(-z)) - 1) / (-z)
    else:
        return 1 / 2

def S(z):
    '''
    Stumpff Function
    '''

    if z > 0:
        return (np.sqrt(z) - np.sin(np.sqrt(z))) / (np.sqrt(z)) ** 3
    elif z < 0:
        return (np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / (np.sqrt(-z)) ** 3
    else:
        return 1 / 6

def newtonRaphson(y, dydx, init_guess, tol, maxiter=1000):
    '''
    Calculate root of single variable function using Newton-Raphson method
    '''
    
    x = init_guess

    for i in range(maxiter):
        x_new = x - y(x) / dydx(x)
        if np.abs(x_new - x) < tol:
            return x_new
        x = x_new
    raise RuntimeError("newtonRaphson: This function did not converge to a solution.")
    
def norm(vec):
    '''
    Vector norm
    '''
    return np.linalg.norm(vec)

def acos(x):
    '''
    Inverse Cosine
    '''
    return np.arccos(np.clip(x, -1.0, 1.0))  # Ensure the value is within the valid range for acos