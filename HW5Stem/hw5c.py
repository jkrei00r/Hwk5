#hw5c.py
#HVS ODEs
# region imports
import numpy as np
from scipy.integrate import solve_ivp  # Import solve_ivp from scipy.integrate
import matplotlib.pyplot as plt
# endregion

# region functions
def ode_system(t, X, *params):
    '''
    The ode system is defined in terms of state variables.
    I have as unknowns:
    x: position of the piston (This is not strictly needed unless I want to know x(t))
    xdot: velocity of the piston
    p1: pressure on right of piston
    p2: pressure on left of the piston
    For initial conditions, we see: x=x0=0, xdot=0, p1=p1_0=p_a, p2=p2_0=p_a
    :param X: The list of state variables.
    :param t: The time for this instance of the function.
    :param params: the list of physical constants for the system.
    :return: The list of derivatives of the state variables.
    '''
    # Unpack the parameters
    A, Cd, ps, pa, V, beta, rho, Kvalve, m, y = params

    # State variables
    x = X[0]  # Position of the piston
    xdot = X[1]  # Velocity of the piston
    p1 = X[2]  # Pressure on the right of the piston
    p2 = X[3]  # Pressure on the left of the piston

    # Calculate derivatives
    xddot = (p1 - p2) * A / m  # Acceleration based on pressure difference
    p1dot = (y * Kvalve * (ps - p1) - rho * A * xdot) * beta / (V * rho)  # Pressure derivative for p1
    p2dot = -(y * Kvalve * (p2 - pa) - rho * A * xdot) * beta / (V * rho)  # Pressure derivative for p2

    # Return the list of derivatives of the state variables
    return [xdot, xddot, p1dot, p2dot]

def main():
    # Define the time array, from 0 to 0.02 seconds, with 200 points
    t = np.linspace(0, 0.02, 200)

    # Parameters for the system
    myargs = (4.909E-4, 0.6, 1.4E7, 1.0E5, 1.473E-4, 2.0E9, 850.0, 2.0E-5, 30, 0.002)

    # Initial conditions: x = 0, xdot = 0, p1 = pa, p2 = pa
    pa = myargs[3]  # pa is the ambient pressure
    ic = [0, 0, pa, pa]  # Initial conditions

    # Call solve_ivp to solve the system of ODEs
    sln = solve_ivp(ode_system, [0, 0.02], ic, args=myargs, t_eval=t)

    # Unpack result into meaningful names
    xvals = sln.y[0]  # Position of the piston
    xdot = sln.y[1]  # Velocity of the piston
    p1 = sln.y[2]  # Pressure on the right of the piston
    p2 = sln.y[3]  # Pressure on the left of the piston

    # Plot the velocity (xdot) as a function of time
    plt.subplot(2, 1, 1)
    plt.plot(t, xdot, 'r-', label=r'$\dot{x}$')  # Use raw string for LaTeX formatting
    plt.title('Velocity as a Function of Time')
    plt.ylabel(r'$\dot{x}$')  # Use raw string for LaTeX formatting
    plt.legend(loc='upper left')

    # Plot p1 and p2 as functions of time
    plt.subplot(2, 1, 2)
    plt.plot(t, p1, 'b-', label=r'$P_1$')  # Use raw string for LaTeX formatting
    plt.plot(t, p2, 'r-', label=r'$P_2$')  # Use raw string for LaTeX formatting
    plt.title('Pressures as a Function of Time')
    plt.xlabel('Time (s)')
    plt.ylabel(r'$P_1, P_2$ (Pa)')  # Use raw string for LaTeX formatting
    plt.legend(loc='upper right')

    # Show the plots
    plt.tight_layout()
    plt.show()

# endregion

# region function calls
if __name__ == "__main__":
    main()
# endregion


