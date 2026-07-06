#sim.py - V0.1

import matplotlib.pyplot as plt 

#Constants
gravity = 9.81 # m/s^2
airDensity = 1.225 # kg/m^3
dragCoefficient = 0.5 # dimensionless
timeStep = 0.1 # seconds

#Inputs
dryMass = 0.1 # kg
radius = 0.001 # m

timeCurve = [0,0.031,0.092,0.139,0.192,0.209,0.231,0.248,0.292,0.37,0.475,0.671,0.702,0.723,0.85,1.063,1.211,1.242,1.303,1.468,1.656,1.821,1.834,1.847,1.86]
thrustCurve = [0,0.946,4.826,9.936,14.09,11.446,7.381,6.151,5.489,4.921,4.448,4.258,4.542,4.164,4.448,4.353,4.353,4.069,4.258,4.353,4.448,4.448,2.933,1.325,0]
massCurve = [10.8,10.782,10.5664,10.1415,9.36163,9.09576,8.84209,8.7012,8.38754,7.89025,7.28777,6.24272,6.07565,5.96368,5.29384,4.14576,3.35665,3.19675,2.88566,2.0155,1.00217,0.103215,0.04445,0.0105492,0]


#Init Variables
time = 0
mass = dryMass + massCurve[0]
altitude = 0
velocity = 0
acceleration = 0
netForce = 0
weight = 0
thrust = 0
drag = 0
avgMassFlowRate = massCurve[0] / timeCurve[len(timeCurve)-1]  # Average mass flow rate over the burn time    
avgThrust = sum(thrustCurve) / timeCurve[len(timeCurve)-1]  # Average thrust over the burn time

#Init Lists
timeList = []
altitudeList = []
velocityList = []


#Functions

def calc_netForce(thrust, weight, drag):
    if velocity >= 0:
        return thrust - weight - drag
    else: 
        return thrust - weight + drag

def calc_acceleration(netForce, mass):
    return netForce / mass

def update_Velocity(velocity, acceleration):
    return velocity + acceleration * timeStep

def update_Altitude(altitude, velocity):
    return altitude + velocity * timeStep

def update_Mass(mass):
    if time > timeCurve[len(timeCurve)-1]:
        return mass
    else:
        return mass - avgMassFlowRate * timeStep

def calc_drag(velocity):
    area = 3.14159 * radius ** 2
    return 0.5 * airDensity * dragCoefficient * area * velocity ** 2

def calc_weight(mass):
    return mass * gravity

def run_simulation():
    global time, mass, altitude, velocity, acceleration, netForce, weight, thrust, drag
    while altitude >= 0:
        # Update time
        time += timeStep

        # Update mass
        mass = update_Mass(mass)

        # Calculate weight
        weight = calc_weight(mass)

        # Calculate drag
        drag = calc_drag(velocity)

        # Determine thrust based on the current time
        if time <= timeCurve[len(timeCurve)-1]:
            thrust = avgThrust  # Use average thrust during the burn time
        else:
            thrust = 0

        # Calculate net force
        netForce = calc_netForce(thrust, weight, drag)

        # Calculate acceleration
        acceleration = calc_acceleration(netForce, mass)

        # Update velocity
        velocity = update_Velocity(velocity, acceleration)

        # Update altitude
        altitude = update_Altitude(altitude, velocity)

        # Store values for plotting
        timeList.append(time)
        altitudeList.append(altitude)
        velocityList.append(velocity)

def plot_results():
    plt.figure(figsize=(12, 6))

    # Plot Altitude vs Time
    plt.subplot(1, 2, 1)
    plt.plot(timeList, altitudeList, label='Altitude (m)', color='blue')
    plt.title('Altitude vs Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Altitude (m)')
    plt.grid()
    plt.legend()

    # Plot Velocity vs Time
    plt.subplot(1, 2, 2)
    plt.plot(timeList, velocityList, label='Velocity (m/s)', color='red')
    plt.title('Velocity vs Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Velocity (m/s)')
    plt.grid()
    plt.legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    run_simulation()
    plot_results()