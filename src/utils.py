'''
Utility functions
'''

import csv
import numpy as np

def stateReader(filespec):
    # This function assumes that the loaded Ephmeris Data is in comma separated value (csv) text format.

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


def lambert(R1,R2,dt,mu,tol=1e-8,maxiter=500,trajectory='pro'):
    '''
    mu         - gravitational parameter (km^3/s^2)
    R1, R2     - initial and final position vectors (km)
    r1, r2     - magnitudes of R1 and R2
    dt         - the time of flight from R1 to R2 (a constant) (s)
    V1, V2     - initial and final velocity vectors (km/s)
    cross12    - cross product of R1 into R2
    dtheta     - angle between R1 and R2
    A          - a constant given by 'Lambert's Solution' lecture notes
    z          - alpha*x^2, where alpha is the reciprocal of the
                semimajor axis and x is the universal anomaly
    y(z)       - a function of z given by Equation 5.38
    F(z,dt)    - a function of the variable z and constant dt,
                - given in 'Lambert's Solution' lecture notes
    dFdz(z)    - the derivative of F(z,t)
    tol        - tolerance on precision of convergence
    f, g       - Lagrange coefficients
    fdot, gdot - time derivative of f and g
    C(z), S(z) - Stumpff functions
    '''

    # Magnitudes of R1 and R2
    r1 = np.linalg.norm(R1)
    r2 = np.linalg.norm(R2)

    # Compute the cross product between R1 and R2
    cross12 = np.cross(R1,R2)

    # Compute the change in true anomaly
    dtheta = np.acos(np.dot(R1,R2)/(r1*r2))

    if trajectory == 'pro':             # For prograde trajectory
        if cross12[2] < 0:
            dtheta = 2*np.pi - dtheta
    elif trajectory == 'retro':         # For retrograde trajectory
        if cross12[2] >= 0:
            dtheta = 2*np.pi - dtheta
    else:
        print(f'Please indicate whether trajectory is prograde or retrograde with keyword argument.')

    # Compute the auxiliary function, A(r1,r2,dtheta)
    A = np.sin(dtheta)*np.sqrt(r1*r2/(1 - np.cos(dtheta)))

    ## Helper functions:
    
    y = lambda z: r1 + r2 + A * (z * S(z) - 1) / np.sqrt(C(z))

    F = lambda z: (y(z) / C(z)) ** 1.5 * S(z) + A * np.sqrt(y(z)) - np.sqrt(mu) * dt

    dFdz = lambda z: (
        np.sqrt(2) / 40 * y(0) ** 1.5 + A / 8 * (np.sqrt(y(0)) + A * np.sqrt(1 / 2 / y(0)))
        if z == 0
        else (y(z) / C(z)) ** 1.5 * (
            1 / 2 / z * (C(z) - 3 * S(z) / 2 / C(z)) + 3 * S(z) ** 2 / 4 / C(z)
        ) + A / 8 * (
            3 * S(z) / C(z) * np.sqrt(y(z)) + A * np.sqrt(C(z) / y(z))
        )
    )

    # Ensure tolerance and maximum iterations are of the correct type
    if not isinstance(tol, (float, int)) or tol <= 0:
        raise ValueError("Tolerance 'tol' must be a positive number.")
    
    # Start with z = 0
    z = 0 

    # Find z by iterating using Newton's method until convergence within the error tolerance:
    z = newtonRaphson(F, dFdz, 0, tol, maxiter)

    # Check the solution
    if np.isnan(z) or np.isinf(z):
        solved = False
    else:
        solved = True 
    
    if solved:
        # Compute the Lagrangian coefficients
        f = 1 - y(z)/r1
        fdot = (np.sqrt(mu)/(r1*r2))*np.sqrt(y(z)/C(z))*(z*S(z)-1)
        g = A*np.sqrt(y(z)/mu)
        gdot = 1 - y(z)/r2

        # Compute the velocities V1 & V2
        V1 = 1/g*(R2 - f*R1)
        V2 = 1/g*(gdot*R2 - R1)
            
        return V1, V2
    else: 
        print('Lambert solver did not converge.')
        return np.array([0, 0, 0]), np.array([0, 0, 0])


def C(z):
    '''
    Stumpff function
    '''
    if z > 0:
        return (1 - np.cos(np.sqrt(z))) / z
    elif z < 0:
        return (np.cosh(np.sqrt(-z)) - 1) / (-z)
    else:
        return 1/2
    
def S(z):
    '''
    Stumpff function
    '''
    if z > 0:
        return (np.sqrt(z) - np.sin(np.sqrt(z))) / (np.sqrt(z))**3
    elif z < 0:
        return (np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / (np.sqrt(-z))**3
    else:
        return 1/6

def newtonRaphson(y,dy,init_guess,tol,maxiter=1000,return_iterations=False):
    '''
    Calculate root of single variable function
	using explicit derivative function
    '''
    
    converged = False # Convergence condition
    x = [init_guess]  # set initial guess
    i = 0             # initialize counter

    # Ensure maximum iterations is of the correct type
    if not isinstance(maxiter, int) or maxiter <= 0:
        raise ValueError("Maximum iterations 'maxiter' must be a positive integer.")

    # begin loop
    while i < maxiter and not(converged):

        # compute new root estimate
        x.append(x[i] - y(x[i]) / dy(x[i]))

        # compute the error
        err = np.abs( (x[i+1] - x[i]) / x[i+1] )

        if err < tol:
            converged = True

        # increment counter
        i = i + 1

    # check for convergence
    if not(converged):
        raise RuntimeError(
            "newtonRaphson: This function did not converge to a solution."
        )

    sol  = x[-1] # final root estimate
    iter_count = i     # total ierations

    if return_iterations:
        return sol, iter_count
    else:
        return sol
    

def norm(vec):
    return np.linalg.norm(vec)