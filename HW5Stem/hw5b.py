#hw5b.py
# Beyond BYOMD: Head Loss Calculator and Moody Diagram Plotter
# region imports
import hw5a as pta  # Import the Moody diagram plotting module (which contains plotMoody)
import random as rnd  # Import random for probabilistic friction factor in transition region
from matplotlib import pyplot as plt  # Import pyplot for additional plotting functionality
import numpy as np  # Import numpy for numerical calculations
# endregion

# region functions

def ffPoint(Re, rr):
    """
    Calculate the friction factor based on the Reynolds number (Re) and relative roughness (rr).
    The function handles laminar, turbulent, and transitional flow regimes.

    :param Re: Reynolds number, dimensionless
    :param rr: Relative roughness (ε/D), dimensionless
    :return: Friction factor (f), dimensionless
    """
    if Re >= 4000:
        # Use Colebrook equation for turbulent flow
        return pta.ff(Re, rr, CBEQN=True)
    elif Re <= 2000:
        # Use laminar flow equation for laminar flow
        return pta.ff(Re, rr)
    else:
        # Transitional flow: Use a probabilistic approach
        CBff = pta.ff(4000, rr, CBEQN=True)  # Friction factor at Re=4000 using Colebrook
        Lamff = pta.ff(2000, rr)  # Friction factor at Re=2000 using laminar equation
        mean = Lamff + (CBff - Lamff) * (Re - 2000) / 2000  # Linear interpolation for mean
        sig = 0.2 * mean  # Standard deviation is 20% of the mean
        return rnd.normalvariate(mean, sig)  # Randomly select from a normal distribution

def PlotPoint(Re, f, plot_points):
    """
    Plot a point on the Moody diagram.

    :param Re: Reynolds number, dimensionless
    :param f: Friction factor, dimensionless
    :param plot_points: List of points to plot
    :return: None (updates the Moody diagram plot)
    """
    plot_points.append((Re, f))  # Append the new point to the list
    pta.plotMoody(plotPoint=True, pt=(Re, f))  # Plot the current point on the Moody diagram
    plt.draw()  # Draw the updated plot
    plt.pause(0.1)  # Pause to allow the plot to update
    return plot_points  # Return the updated list of points

def calculate_head_loss(diameter_inches, roughness_micro_inches, flow_rate_gpm):
    """
    Calculate the head loss per foot (hf/L) in English units based on user inputs.

    :param diameter_inches: Pipe diameter in inches
    :param roughness_micro_inches: Pipe roughness in micro-inches (10^-6 inches)
    :param flow_rate_gpm: Flow rate in gallons per minute (gpm)
    :return: Head loss per foot (hf/L) in feet of fluid per foot of pipe
    """
    # Convert inputs to consistent units
    diameter_feet = diameter_inches / 12  # Convert diameter to feet
    roughness_feet = roughness_micro_inches * 1e-6 / 12  # Convert roughness to feet
    flow_rate_cfs = flow_rate_gpm * 0.002228  # Convert flow rate to cubic feet per second (cfs)

    # Calculate Reynolds number (Re)
    kinematic_viscosity = 1.08e-5  # Kinematic viscosity of water at 60°F in ft^2/s
    velocity = flow_rate_cfs / (np.pi * (diameter_feet / 2)**2)  # Velocity in ft/s
    Re = velocity * diameter_feet / kinematic_viscosity

    # Calculate relative roughness (ε/D)
    rr = roughness_feet / diameter_feet

    # Calculate friction factor (f)
    f = ffPoint(Re, rr)

    # Calculate head loss per foot (hf/L) using Darcy-Weisbach equation
    g = 32.2  # Acceleration due to gravity in ft/s^2
    hf_L = f * (velocity**2) / (2 * g * diameter_feet)

    return hf_L, Re, f

def main():
    """
    Main function to interact with the user, calculate head loss, and plot results on the Moody diagram.
    """
    # Create a list to store plot points
    plot_points = []

    # Initialize the plot
    plt.ion()  # Turn on interactive mode
    pta.plotMoody()  # Draw the initial empty Moody diagram

    while True:
        # Prompt the user for input data
        diameter_inches = float(input("Enter the pipe diameter in inches: "))
        roughness_micro_inches = float(input("Enter the pipe roughness in micro-inches: "))
        flow_rate_gpm = float(input("Enter the flow rate in gallons per minute: "))

        # Calculate head loss and friction factor
        hf_L, Re, f = calculate_head_loss(diameter_inches, roughness_micro_inches, flow_rate_gpm)

        # Determine flow type
        if Re < 2000:
            flow_type = 'laminar'
        elif Re > 4000:
            flow_type = 'turbulent'
        else:
            flow_type = 'transition'

        # Display results
        print(f"Head loss per foot (hf/L): {hf_L:.6f} ft/ft")
        print(f"Reynolds number (Re): {Re:.2f}")
        print(f"Friction factor (f): {f:.4f}")
        print(f"Flow type: {flow_type}")

        # Plot the point on the Moody diagram and retain previous points
        plot_points = PlotPoint(Re, f, plot_points)

        # Ask the user if they want to continue entering more parameters
        continue_input = input("Do you want to enter another set of parameters? (yes/no): ").lower()
        if continue_input != 'yes':
            break

    # Final plot display
    plt.ioff()  # Turn off interactive mode
    plt.show()  # Display the final plot

# endregion

# region function calls
if __name__ == "__main__":
    main()
# endregion
