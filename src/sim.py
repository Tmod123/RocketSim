#sim.py - V1.6

from ast import While

import matplotlib.pyplot as plt 
import numpy as np
import csv

def import_Motor_Data(filename):
    times = []
    thrusts = []
    times.append(0.000)
    thrusts.append(0.000)
    with open(filename, newline="") as f:
        reader = csv.reader(f)

        for row in reader:
            if not row:
                continue

            try:
                times.append(float(row[0]))
                thrusts.append(float(row[1]))
            except (ValueError, IndexError):
                print(f"Invalid value in row: {row}")

        return [times, thrusts]

#print(import_Motor_Data("Estes_E12.csv"))

#Constants
gravity = 9.81 # m/s^2
dragCoefficient = 0.5 # dimensionless
timeStep = 0.05 # seconds

#Inputs
dryMass = 0.0605 # kg
fuelMass = 0.011 # kg
radius = 0.0127 # m
parachuteArea = 0.0707 # m^2
parachuteDragCoefficient = 0.80 # dimensionless
rodLength = 1.0 # m
#fileName = input("Enter the CSV file containing motor data: ") # CSV file containing motor data

#motor_csv = import_Motor_Data(fileName)
motor_csv = import_Motor_Data("Motors/Estes_C6.csv")

timeCurve = motor_csv[0]
thrustCurve = motor_csv[1]
#timeCurve = [0,0.031,0.092,0.139,0.192,0.209,0.231,0.248,0.292,0.37,0.475,0.671,0.702,0.723,0.85,1.063,1.211,1.242,1.303,1.468,1.656,1.821,1.834,1.847,1.86]
#thrustCurve = [0,0.946,4.826,9.936,14.09,11.446,7.381,6.151,5.489,4.921,4.448,4.258,4.542,4.164,4.448,4.353,4.353,4.069,4.258,4.353,4.448,4.448,2.933,1.325,0]
#massCurve = [10.8,10.782,10.5664,10.1415,9.36163,9.09576,8.84209,8.7012,8.38754,7.89025,7.28777,6.24272,6.07565,5.96368,5.29384,4.14576,3.35665,3.19675,2.88566,2.0155,1.00217,0.103215,0.04445,0.0105492,0]


#Init Variables
time = 0
mass = dryMass + (fuelMass)  # Initial mass in kg (dry mass + initial propellant mass)
altitude = 0
velocity = 0
acceleration = 0
netForce = 0
weight = 0
thrust = 0
drag = 0
#avgMassFlowRate = (massCurve[0] / timeCurve[len(timeCurve)-1])/1000  # Average mass flow rate over the burn time    

#Init Lists
timeList = []
massList = []
altitudeList = []
velocityList = []
accelerationList = []
netForceList = []
weightList = []
thrustList = []
dragList = []
densityList = []



#Functions

def acceleration(time, altitude, velocity):
    return netForce(time, velocity, altitude) / mass(time)

def netForce(time, velocity, altitude):
    if velocity >= 0:
        return thrust(time) - weight(time) - drag(time, velocity, altitude)
    else:
        return thrust(time) - weight(time) + drag(time, velocity, altitude)

def thrust(time):
    if time > timeCurve[-1]:
        return 0
    elif time <= timeCurve[0]:
        return thrustCurve[1]*time/timeCurve[1]
    else:
        for i in range(len(timeCurve)-1):
            if timeCurve[i] < time <= timeCurve[i+1]:
                Time1= timeCurve[i]
                Time2= timeCurve[i+1]
                Thrust1= thrustCurve[i]
                Thrust2= thrustCurve[i+1]
                break
        return Thrust1 + (Thrust2 - Thrust1) * (time - Time1) / (Time2 - Time1)


def weight(time):
    return mass(time) * gravity

def drag(time, velocity, altitude):
    if velocity >= 0:
        area = 3.14159 * radius ** 2
        return 0.5 * air_density(altitude) * dragCoefficient * area * velocity ** 2
    else:
        return 0.5 * air_density(altitude) * parachuteDragCoefficient * parachuteArea * velocity ** 2

#def mass(time):
#    if time > timeCurve[len(timeCurve)-1]:
#        return dryMass
#    else:
#        return dryMass + (massCurve[0]/1000) - avgMassFlowRate * time

def mass(time):
    currentImpulse=0
    netImpulse=0
    for i in range(1,len(timeCurve)):
        if time < timeCurve[0]:
            currentImpulse = thrust(time)*time
        elif timeCurve[i] <= time:
            currentImpulse += (thrustCurve[i]*(timeCurve[i]-timeCurve[i-1]))
        else:
            currentImpulse += (thrust(time)*(time-timeCurve[i-1]))
            break
    for i in range(1,len(timeCurve)):
        netImpulse += (thrustCurve[i]*(timeCurve[i]-timeCurve[i-1]))
    currentFuelMass=fuelMass*(1-(currentImpulse/netImpulse))
    return dryMass + currentFuelMass


def air_density(altitude):
    return 1.225 * (1-((0.0065/288.15)*altitude))**((gravity/(287.05*0.0065))-1)
    
def rk4_step(time, altitude, velocity):
    k1_a = velocity
    k1_v = acceleration(time, altitude, velocity)

    k2_a = velocity + k1_v * timeStep/2
    k2_v = acceleration(time + timeStep/2, altitude + k1_a * timeStep/2, velocity + k1_v * timeStep/2)

    k3_a = velocity + k2_v * timeStep/2
    k3_v = acceleration(time + timeStep/2, altitude + k2_a * timeStep/2, velocity + k2_v * timeStep/2)

    k4_a = velocity + k3_v * timeStep
    k4_v = acceleration(time + timeStep, altitude + k3_a * timeStep, velocity + k3_v * timeStep)

    new_altitude = altitude + (timeStep/6) * (k1_a + 2*k2_a + 2*k3_a + k4_a)
    new_velocity = velocity + (timeStep/6) * (k1_v + 2*k2_v + 2*k3_v + k4_v)

    return new_altitude, new_velocity
    
def run_simulation():
    global time, mass, altitude, velocity, acceleration, netForce, weight, thrust, drag
    while altitude >= -0.01:
        timeList.append(time)
        altitudeList.append(altitude)
        velocityList.append(velocity)
        massList.append(mass(time))
        accelerationList.append(acceleration(time, altitude, velocity))
        netForceList.append(netForce(time, velocity, altitude))
        weightList.append(weight(time))
        thrustList.append(thrust(time))
        dragList.append(drag(time, velocity, altitude))
        densityList.append(air_density(altitude))
        new_altitude, new_velocity = rk4_step(time, altitude, velocity)
        altitude = new_altitude
        velocity = new_velocity
        time += timeStep
        

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
   for i in range(len(timeList)):
       print(f"Time: {timeList[i]:.2f} s, Mass: {massList[i]:.2f} kg, Altitude: {altitudeList[i]:.2f} m, Velocity: {velocityList[i]:.2f} m/s, acceleration: {accelerationList[i]:.2f} m/s^2, Net Force: {netForceList[i]:.2f} N, Weight: {weightList[i]:.2f} N, Thrust: {thrustList[i]:.2f} N, Drag: {dragList[i]:.2f} N, Air Density: {densityList[i]:.2f} kg/m^3")
   for i in range(len(altitudeList)-1):
       if altitudeList[i] < rodLength <= altitudeList[i+1]:
           rodVelocity = np.sqrt(velocityList[i]**2 + 2 * accelerationList[i] * (rodLength - altitudeList[i]))
           break
   print(f"Velocity off rod: {rodVelocity:.2f} m/s, Apogee: {max(altitudeList):.2f} m, Max Velocity: {max(velocityList):.2f} m/s, Max acceleration: {max(accelerationList):.2f} m/s^2, Time to Apogee: {timeList[altitudeList.index(max(altitudeList))]:.2f} s, Flight Time: {timeList[-1]:.2f} s, Ground hit velocity: {velocityList[-1]:.2f} m/s")
   plot_results()

#print(mass(0.5))