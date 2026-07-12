#sim.py - V2.01

import matplotlib.pyplot as plt 
import numpy as np
import csv

def import_Motor_Data(filename):
    times = [0.000]
    thrusts = [0.000]
    with open(filename, newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            try:
                times.append(float(row[0]))
                thrusts.append(float(row[1]))
            except (ValueError, IndexError):
                pass
        return [times, thrusts]

#Constants
gravity = 9.81 # m/s^2
timeStep = 0.05 # seconds

#Inputs
dryMass = 0.0605 # kg
fuelMass = 0.011 # kg
radius = 0.0127 # m
dragCoefficient = 0.5 # dimensionless
parachuteArea = 0.0707 # m^2
parachuteDragCoefficient = 0.80 # dimensionless
rodLength = 1.0 # m
launchAngle = np.radians(5) #radians of x degrees
windSpeed = 0 # m/s
#fileName = input("Enter the CSV file containing motor data: ") # CSV file containing motor data

#motor_csv = import_Motor_Data(fileName)
motor_csv = import_Motor_Data("Motors/Estes_C6.csv")
timeCurve = motor_csv[0]
thrustCurve = motor_csv[1]

#Init Variables
time = 0
y_pos = 0
y_velocity = 0
x_pos = 0
x_velocity = 0
leftRod = False

#Init Lists
timeList = []
massList = []
x_posList = []
y_posList = []
x_velocityList = []
y_velocityList = []
x_accelerationList = []
y_accelerationList = []
x_netForceList = []
y_netForceList = []
weightList = []
x_thrustList = []
y_thrustList = []
x_dragList = []
y_dragList = []
densityList = []



#Functions

def air_density(y_pos):
    return 1.225 * (1-((0.0065/288.15)*y_pos))**((gravity/(287.05*0.0065))-1)

def netSpeed(x_velocity, y_velocity):
    return np.sqrt((x_velocity**2) + (y_velocity**2))

def on_rod(x_pos, y_pos):
    global leftRod
    if leftRod:
        return False
    if np.sqrt((x_pos**2) + (y_pos**2)) <= rodLength:
        return True
    leftRod = True
    return False

def flight_angle(x_pos, y_pos, x_velocity, y_velocity):
    if(on_rod(x_pos, y_pos)):
        return launchAngle
    return np.arctan2(x_velocity, y_velocity)

def thrust_magnitude(time):
    if time > timeCurve[-1]:
        return 0
    if time <= timeCurve[0]:
        return (thrustCurve[1]*time/timeCurve[1])
    for i in range(len(timeCurve)-1):
        if timeCurve[i] < time <= timeCurve[i+1]:
            Time1= timeCurve[i]
            Time2= timeCurve[i+1]
            Thrust1= thrustCurve[i]
            Thrust2= thrustCurve[i+1]
            return (Thrust1 + (Thrust2 - Thrust1) * (time - Time1) / (Time2 - Time1))
    return 0

def x_component(magnitude, x_pos, y_pos, x_velocity, y_velocity):
    return magnitude * np.sin(flight_angle(x_pos, y_pos, x_velocity, y_velocity))

def y_component(magnitude, x_pos, y_pos, x_velocity, y_velocity):
    return magnitude * np.cos(flight_angle(x_pos, y_pos, x_velocity, y_velocity))

def acceleration_x(time, x_pos, y_pos, x_velocity, y_velocity):
    return netForce_x(time, x_pos, y_pos, x_velocity, y_velocity) / mass(time)

def acceleration_y(time, x_pos, y_pos, x_velocity, y_velocity):
    return netForce_y(time, x_pos, y_pos, x_velocity, y_velocity) / mass(time)

def netForce_x(time, x_pos, y_pos, x_velocity, y_velocity):
    if y_pos == 0:
        return x_component(thrust_magnitude(time), x_pos, y_pos, x_velocity, y_velocity)
    else:
        return x_component(thrust_magnitude(time), x_pos, y_pos, x_velocity, y_velocity) - x_component(drag_magnitude(x_pos, y_pos, x_velocity, y_velocity),x_pos, y_pos, x_velocity, y_velocity)    

def netForce_y(time, x_pos, y_pos, x_velocity, y_velocity):
    if y_pos == 0:
        return y_component(thrust_magnitude(time), x_pos, y_pos, x_velocity, y_velocity)
    else:
        return y_component(thrust_magnitude(time), x_pos, y_pos, x_velocity, y_velocity) - weight(time) - y_component(drag_magnitude(x_pos, y_pos, x_velocity, y_velocity),x_pos, y_pos, x_velocity, y_velocity)

def drag_magnitude(x_pos, y_pos, x_velocity, y_velocity):
    speed = netSpeed(x_velocity, y_velocity)
    if y_velocity >= 0:
        A, cd = 3.14159 * radius ** 2, dragCoefficient        
    else:
        A, cd = parachuteArea, parachuteDragCoefficient
    return (0.5 * air_density(y_pos) * cd * A * speed ** 2)

def weight(time):
    return mass(time) * gravity

def mass(time):
    currentImpulse=0
    netImpulse=0
    for i in range(1,len(timeCurve)):
        if time < timeCurve[0]:
            currentImpulse = thrust_magnitude(time) * time
        elif timeCurve[i] <= time:
            currentImpulse += 0.5*(thrustCurve[i]+thrustCurve[i-1])*(timeCurve[i]-timeCurve[i-1])
        else:
            currentImpulse += 0.5*(thrust_magnitude(time)+thrustCurve[i-1]) * (time-timeCurve[i-1])
            break
    for i in range(1,len(timeCurve)):
        netImpulse += (0.5*(thrustCurve[i]+thrustCurve[i-1])*(timeCurve[i]-timeCurve[i-1]))
    currentFuelMass=fuelMass*(1-(currentImpulse/netImpulse))
    return dryMass + currentFuelMass

def rk4_step(time, x_pos, y_pos, x_velocity, y_velocity):
    k1_x = x_velocity
    k1_y = y_velocity
    k1_vx = acceleration_x(time, x_pos, y_pos, x_velocity, y_velocity)
    k1_vy = acceleration_y(time, x_pos, y_pos, x_velocity, y_velocity)

    k2_x = x_velocity + k1_vx * timeStep/2
    k2_y = y_velocity + k1_vy * timeStep/2
    k2_vx = acceleration_x(time + timeStep/2, x_pos + k1_x * timeStep/2, y_pos + k1_y * timeStep/2, x_velocity + k1_vx * timeStep/2, y_velocity + k1_vy * timeStep/2)
    k2_vy = acceleration_y(time + timeStep/2, x_pos + k1_x * timeStep/2, y_pos + k1_y * timeStep/2, x_velocity + k1_vx * timeStep/2, y_velocity + k1_vy * timeStep/2)

    k3_x = x_velocity + k2_vx * timeStep/2
    k3_y = y_velocity + k2_vy * timeStep/2
    k3_vx = acceleration_x(time + timeStep/2, x_pos + k2_x * timeStep/2, y_pos + k2_y * timeStep/2, x_velocity + k2_vx * timeStep/2, y_velocity + k2_vy * timeStep/2)
    k3_vy = acceleration_y(time + timeStep/2, x_pos + k2_x * timeStep/2, y_pos + k2_y * timeStep/2, x_velocity + k2_vx * timeStep/2, y_velocity + k2_vy * timeStep/2)

    k4_x = x_velocity + k3_vx * timeStep
    k4_y = y_velocity + k3_vy * timeStep
    k4_vx = acceleration_x(time + timeStep, x_pos + k3_x * timeStep, y_pos + k3_y * timeStep, x_velocity + k3_vx * timeStep, y_velocity + k3_vy * timeStep)
    k4_vy = acceleration_y(time + timeStep, x_pos + k3_x * timeStep, y_pos + k3_y * timeStep, x_velocity + k3_vx * timeStep, y_velocity + k3_vy * timeStep)

    new_x_pos = x_pos + (timeStep/6) * (k1_x + 2*k2_x + 2*k3_x + k4_x)
    new_y_pos = y_pos + (timeStep/6) * (k1_y + 2*k2_y + 2*k3_y + k4_y)
    new_x_velocity = x_velocity + (timeStep/6) * (k1_vx + 2*k2_vx + 2*k3_vx + k4_vx)
    new_y_velocity = y_velocity + (timeStep/6) * (k1_vy + 2*k2_vy + 2*k3_vy + k4_vy)

    return new_x_pos, new_y_pos, new_x_velocity, new_y_velocity
    
def run_simulation():
    global time, x_pos, y_pos, x_velocity, y_velocity
    while y_pos >= -0.01:

        current_mass = mass(time)
        current_thrust_x = x_component(thrust_magnitude(time), x_pos, y_pos, x_velocity, y_velocity)
        current_thrust_y = y_component(thrust_magnitude(time), x_pos, y_pos, x_velocity, y_velocity)
        current_weight = current_mass * gravity
        current_drag_x = x_component(drag_magnitude(x_pos, y_pos, x_velocity, y_velocity), x_pos, y_pos, x_velocity, y_velocity)
        current_drag_y = y_component(drag_magnitude(x_pos, y_pos, x_velocity, y_velocity), x_pos, y_pos, x_velocity, y_velocity)

        if y_pos == 0:
            current_netForce_x = current_thrust_x
            current_netForce_y = current_thrust_y 
        else:
            current_netForce_x = current_thrust_x - current_drag_x
            current_netForce_y = current_thrust_y - current_weight - current_drag_y

        current_acceleration_x = current_netForce_x / current_mass
        current_acceleration_y = current_netForce_y / current_mass

        timeList.append(time)
        x_posList.append(x_pos)
        y_posList.append(y_pos)
        x_velocityList.append(x_velocity)  
        y_velocityList.append(y_velocity)
        massList.append(current_mass)
        x_accelerationList.append(current_acceleration_x)
        y_accelerationList.append(current_acceleration_y)
        x_netForceList.append(current_netForce_x)
        y_netForceList.append(current_netForce_y)
        weightList.append(current_weight)
        x_thrustList.append(current_thrust_x)
        y_thrustList.append(current_thrust_y)
        x_dragList.append(current_drag_x)
        y_dragList.append(current_drag_y)
        densityList.append(air_density(y_pos))
        new_x_pos, new_y_pos, new_x_velocity, new_y_velocity = rk4_step(time, x_pos, y_pos, x_velocity, y_velocity)
        x_pos = new_x_pos
        y_pos = new_y_pos
        x_velocity = new_x_velocity
        y_velocity = new_y_velocity
        time += timeStep
        

def plot_results():
    fig, ax1 = plt.subplots(figsize=(10, 8))
    ax1.plot(timeList, y_posList, color='blue', label='Altitude (m)')
    ax1.set_ylabel('Altitude (m)', color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')

    ax2 = ax1.twinx()

    ax2.plot(timeList, y_velocityList, color='red', label='Velocity (m/s)')
    ax2.set_ylabel('Velocity (m/s)', color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

    plt.title('Vertical Motion vs Time')
    plt.tight_layout()
    plt.show()
    

if __name__ == "__main__":
    run_simulation()
    for i in range(len(timeList)):
        print(f"Time: {timeList[i]:.2f} s, Mass: {massList[i]:.2f} kg, Altitude: {y_posList[i]:.2f} m, Velocity: {y_velocityList[i]:.2f} m/s, acceleration: {y_accelerationList[i]:.2f} m/s^2, Net Force: {y_netForceList[i]:.2f} N, Weight: {weightList[i]:.2f} N, Thrust: {y_thrustList[i]:.2f} N, Drag: {y_dragList[i]:.2f} N, Air Density: {densityList[i]:.2f} kg/m^3")
    rodVelocity = 0
    for i in range(len(y_posList)-1):
        dist1=np.sqrt(y_posList[i]**2+x_posList[i]**2)
        dist2=np.sqrt(y_posList[i+1]**2+x_posList[i+1]**2)
        if dist1 < rodLength <= dist2:
            rodVelocity = np.sqrt((y_velocityList[i]**2 + 2 * y_accelerationList[i] * (rodLength*np.cos(launchAngle) - y_posList[i]))+(x_velocityList[i]**2 + 2 * x_accelerationList[i] * (rodLength*np.sin(launchAngle) - x_posList[i])))
            break
    print(f"Velocity off rod: {rodVelocity:.2f} m/s, Apogee: {max(y_posList):.2f} m, Max Velocity: {max(y_velocityList):.2f} m/s, Max acceleration: {max(y_accelerationList):.2f} m/s^2, Time to Apogee: {timeList[y_posList.index(max(y_posList))]:.2f} s, Flight Time: {timeList[-1]:.2f} s, Ground hit velocity: {y_velocityList[-1]:.2f} m/s, Range from launch: {x_posList[-1]:.2f} m")
    plot_results()
