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
    Swamee-Jain approximation for the friction factor in turbulent flow.
    This is used as an initial guess for the Colebrook equation solver.

    :param Re: Reynolds number
    :param rr: relative roughness
    :return: approximate friction factor
    """
    return 0.25 / (np.log10((rr / 3.7) + (5.74 / Re ** 0.9))) ** 2


def ff(Re, rr, CBEQN=False):
    """
    This function calculates the friction factor for a pipe based on the
    notion of laminar, turbulent, and transitional flow.

    :param Re: the Reynolds number under question.
    :param rr: the relative pipe roughness (expect between 0 and 0.05)
    :param CBEQN: boolean to indicate if Colebrook equation (True) or laminar equation should be used
    :return: the (Darcy) friction factor
    """
    if CBEQN:
        # Colebrook equation for turbulent flow
        def colebrook(f):
            if f <= 0 or Re <= 0:
                return 1e6  # Return a large value to avoid invalid operations
            return 1 / np.sqrt(f) + 2.0 * np.log10((rr / 3.7) + (2.51 / (Re * np.sqrt(f))))

        # Initial guess using Swamee-Jain approximation
        initial_guess = swamee_jain(Re, rr)
        try:
            result = fsolve(colebrook, initial_guess, full_output=True)
            if result[2] == 1:  # Check if fsolve converged
                return result[0][0]
            else:
                # Fallback to Swamee-Jain approximation if fsolve fails
                return swamee_jain(Re, rr)
        except:
            # Fallback to Swamee-Jain approximation if fsolve fails
            return swamee_jain(Re, rr)
    else:
        # Laminar flow equation
        return 64 / Re


def plotMoody(plotPoint=False, pt=(0, 0)):
    """
    This function produces the Moody diagram for a Re range from 1 to 10^8 and
    for relative roughness from 0 to 0.05 (20 steps). The laminar region is described
    by the simple relationship of f=64/Re whereas the turbulent region is described by
    the Colebrook equation.

    :param plotPoint: boolean to indicate if a specific point should be plotted
    :param pt: tuple representing the specific point to plot (Re, f)
    :return: just shows the plot, nothing returned
    """
    # Step 1: Create logspace arrays for ranges of Re
    ReValsCB = np.logspace(np.log10(4000), np.log10(1e8),
                           200)  # for use with Colebrook equation (i.e., Re in range from 4000 to 10^8)
    ReValsL = np.logspace(np.log10(600.0), np.log10(2000.0),
                          20)  # for use with Laminar flow (i.e., Re in range from 600 to 2000)
    ReValsTrans = np.logspace(np.log10(2000.0), np.log10(4000.0),
                              20)  # for use with Transition flow (i.e., Re in range from 2000 to 4000)

    # Step 2: Create array for range of relative roughnesses
    rrVals = np.array(
        [0, 1E-6, 5E-6, 1E-5, 5E-5, 1E-4, 2E-4, 4E-4, 6E-4, 8E-4, 1E-3, 2E-3, 4E-3, 6E-3, 8E-8, 1.5E-2, 2E-2, 3E-2,
         4E-2, 5E-2])

    # Step 3: Calculate the friction factor in the laminar range
    ffLam = np.array([ff(Re, 0) for Re in ReValsL])  # use list comprehension for all Re in ReValsL and calling ff
    ffTrans = np.array(
        [ff(Re, 0) for Re in ReValsTrans])  # use list comprehension for all Re in ReValsTrans and calling ff

    # Step 4: Calculate friction factor values for each rr at each Re for turbulent range
    ffCB = np.array([[ff(Re, rr, True) for Re in ReValsCB] for rr in rrVals])

    # Step 5: Construct the plot
    plt.loglog(ReValsL, ffLam, 'b-', label='Laminar')  # plot the laminar part as a solid line
    plt.loglog(ReValsTrans, ffTrans, 'b--', label='Transition')  # plot the transition part as a dashed line
    for nRelR in range(len(ffCB)):
        plt.loglog(ReValsCB, ffCB[nRelR], color='k',
                   label=f'rr={rrVals[nRelR]}')  # plot the lines for the turbulent region for each pipe roughness
        plt.annotate(xy=(1e8, ffCB[nRelR][-1]), text=f'{rrVals[nRelR]}',
                     fontsize=8)  # put a label at end of each curve on the right

    plt.xlim(600, 1e8)
    plt.ylim(0.008, 0.10)
    plt.xlabel(r"Reynolds number $Re$", fontsize=16)
    plt.ylabel(r"Friction factor $f$", fontsize=16)
    plt.text(2.5e8, 0.02, r"Relative roughness $\frac{\epsilon}{d}$", rotation=90, fontsize=16)
    ax = plt.gca()  # capture the current axes for use in modifying ticks, grids, etc.
    ax.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=12)  # format tick marks
    ax.tick_params(axis='both', grid_linewidth=1, grid_linestyle='solid', grid_alpha=0.5)
    ax.tick_params(axis='y', which='minor')
    ax.yaxis.set_minor_formatter(FormatStrFormatter("%.3f"))
    plt.grid(which='both')
    if plotPoint:
        plt.plot(pt[0], pt[1], 'ro', markersize=12, markeredgecolor='red', markerfacecolor='none')

    plt.show()


def main():
    plotMoody()


# endregion

# region function calls
if __name__ == "__main__":
    main()
# endregion