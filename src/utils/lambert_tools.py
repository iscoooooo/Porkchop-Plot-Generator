'''
Lambert's Problem Tools
'''

# Third-party Libraries
import numpy as np

# Porkchop-Plot-Generator Libraries
from utils.numerical_tools import newtonRaphson


def lambert_solver(R1, R2, dt, mu, tol=1e-6, maxiter=10000, trajectory='pro'):
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