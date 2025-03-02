#hw5a.py
# BYOMD (Build Your Own Moody Diagram)
# region imports
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

# endregion

# region functions

def swamee_jain(Re, rr):
    """
    Calculate the approximate friction factor using the Swamee-Jain equation.
    This is used as an initial guess for the Colebrook equation solver and as a fallback.

    :param Re: Reynolds number, dimensionless
    :param rr: Relative roughness (ε/D), dimensionless
    :return: Approximate friction factor (f), dimensionless
    """
    # Swamee-Jain approximation for turbulent flow
    return 0.25 / (np.log10((rr / 3.7) + (5.74 / Re ** 0.9))) ** 2


def ff(Re, rr, CBEQN=False):
    """
    Calculate the Darcy friction factor for pipe flow based on the Reynolds number and relative roughness.
    The function handles both laminar and turbulent flow regimes.

    :param Re: Reynolds number, dimensionless
    :param rr: Relative roughness (ε/D), dimensionless
    :param CBEQN: Boolean flag to indicate whether to use the Colebrook equation (True) or the laminar flow equation (False)
    :return: Darcy friction factor (f), dimensionless
    """
    if CBEQN:
        # Colebrook equation for turbulent flow
        def colebrook(f):
            """
            Colebrook equation to solve for the friction factor in turbulent flow.
            The equation is implicit and requires numerical methods to solve.

            :param f: Friction factor, dimensionless
            :return: Residual of the Colebrook equation
            """
            if f <= 0 or Re <= 0:
                return 1e6  # Return a large value to avoid invalid operations
            return 1 / np.sqrt(f) + 2.0 * np.log10((rr / 3.7) + (2.51 / (Re * np.sqrt(f))))

        # Initial guess using Swamee-Jain approximation
        initial_guess = swamee_jain(Re, rr)
        try:
            # Use fsolve to numerically solve the Colebrook equation
            result = fsolve(colebrook, initial_guess, full_output=True)
            if result[2] == 1:  # Check if fsolve converged
                return result[0][0]
            else:
                # Fallback to Swamee-Jain approximation if fsolve fails to converge
                return swamee_jain(Re, rr)
        except:
            # Fallback to Swamee-Jain approximation if fsolve fails
            return swamee_jain(Re, rr)
    else:
        # Laminar flow equation (f = 64 / Re)
        return 64 / Re


def plotMoody(plotPoint=False, pt=(0, 0)):
    """
    Generate and display a Moody diagram, which plots the Darcy friction factor (f) against the Reynolds number (Re)
    for various relative roughness values (ε/D). The diagram covers laminar, transitional, and turbulent flow regimes.

    :param plotPoint: Boolean flag to indicate whether to plot a specific point on the diagram
    :param pt: Tuple (Re, f) representing the specific point to plot
    :return: None (displays the plot)
    """
    # Step 1: Create logspace arrays for ranges of Reynolds numbers (Re)
    ReValsCB = np.logspace(np.log10(4000), np.log10(1e8), 200)  # Turbulent flow range (Re from 4000 to 10^8)
    ReValsL = np.logspace(np.log10(600.0), np.log10(2000.0), 20)  # Laminar flow range (Re from 600 to 2000)
    ReValsTrans = np.logspace(np.log10(2000.0), np.log10(4000.0), 20)  # Transitional flow range (Re from 2000 to 4000)

    # Step 2: Create an array for a range of relative roughness values (ε/D)
    rrVals = np.array(
        [0, 1E-6, 5E-6, 1E-5, 5E-5, 1E-4, 2E-4, 4E-4, 6E-4, 8E-4, 1E-3, 2E-3, 4E-3, 6E-3, 8E-8, 1.5E-2, 2E-2, 3E-2,
         4E-2, 5E-2])

    # Step 3: Calculate the friction factor in the laminar range using the laminar flow equation
    ffLam = np.array([ff(Re, 0) for Re in ReValsL])  # Use list comprehension to calculate f for all Re in ReValsL

    # Step 4: Calculate the friction factor in the transitional range
    ffTrans = np.array([ff(Re, 0) for Re in ReValsTrans])  # Use list comprehension to calculate f for all Re in ReValsTrans

    # Step 5: Calculate the friction factor values for each relative roughness (rr) at each Reynolds number (Re) for the turbulent range
    ffCB = np.array([[ff(Re, rr, True) for Re in ReValsCB] for rr in rrVals])  # Use list comprehension to calculate f for each rr

    # Step 6: Create the plot with the calculated values
    plt.loglog(ReValsL, ffLam, 'b-', label='Laminar')  # Plot the laminar region as a solid blue line
    plt.loglog(ReValsTrans, ffTrans, 'b--', label='Transition')  # Plot the transition region as a dashed blue line
    for nRelR in range(len(ffCB)):
        plt.loglog(ReValsCB, ffCB[nRelR], color='k', label=f'rr={rrVals[nRelR]}')  # Plot the turbulent region for each relative roughness
        plt.annotate(xy=(1e8, ffCB[nRelR][-1]), text=f'{rrVals[nRelR]}', fontsize=8)  # Annotate each turbulent curve at the end

    # Step 7: Set the plot limits, labels, and formatting
    plt.xlim(600, 1e8)  # Set x-axis limits for Reynolds number
    plt.ylim(0.008, 0.10)  # Set y-axis limits for friction factor
    plt.xlabel(r"Reynolds number $Re$", fontsize=16)  # Label for x-axis
    plt.ylabel(r"Friction factor $f$", fontsize=16)  # Label for y-axis
    plt.text(2.5e8, 0.02, r"Relative roughness $\frac{\epsilon}{d}$", rotation=90, fontsize=16)  # Label for relative roughness

    # Step 8: Format the plot axes and grid
    ax = plt.gca()  # Get the current axes for customization
    ax.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=12)  # Customize tick marks
    ax.tick_params(axis='both', grid_linewidth=1, grid_linestyle='solid', grid_alpha=0.5)  # Customize grid lines
    ax.tick_params(axis='y', which='minor')  # Minor ticks on y-axis
    ax.yaxis.set_minor_formatter(FormatStrFormatter("%.3f"))  # Format the minor ticks on y-axis
    plt.grid(which='both')  # Enable both major and minor grid lines

    # Step 9: Plot a specific point if requested
    if plotPoint:
        plt.plot(pt[0], pt[1], 'ro', markersize=12, markeredgecolor='red', markerfacecolor='none')  # Plot a red circle for the point
        plt.draw()  # Redraw the plot with the new point
        plt.pause(0.1)  # Pause briefly to update the plot without blocking execution

    # Step 10: Display the plot (this is required to actually show the graph)
    plt.show()  # This will display the plot in a window


def main():
    """
    Main function to generate and display the Moody diagram.
    This function is responsible for generating the plot and displaying it without blocking
    the rest of the program's execution.
    """
    # Generate and display the Moody diagram
    plotMoody()


# endregion

# region function calls
if __name__ == "__main__":
    main()
# endregion
