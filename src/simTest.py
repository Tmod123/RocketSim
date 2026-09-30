#sim.py - V2.3
import matplotlib.pyplot as plt 
import numpy as np
import csv
import sys
import random
import json


def import_Motor_Data(filename):
    times = [0.000]
    thrusts = [0.000]
    with open(filename, newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            if not row:
                continue
            try:
                if float(row[1]) != 0.000:
                    times.append(float(row[0]))
                    thrusts.append(float(row[1]))
            except (ValueError, IndexError):
                pass
        return [times, thrusts]



dryMass = 3.693  # kg — mass properties computed during simulation
fuelMass = 1.552   # kg — motor propellant mass
dragCoefficient = 0.566  # dimensionless — simulation/CFD result
rodLength   = 1.0                  # m — launch condition, not part of <rocket>
launchAngle = np.radians(5)
launchDirection = np.radians(45) #heading
launch_vector = np.array([np.sin(launchAngle)*np.cos(launchDirection), np.sin(launchAngle) * np.sin(launchDirection), np.cos(launchAngle)])
CG_dry = 0.9964  # m from nose tip — computed mass property
CG_wet = 1.0237  # m from nose tip — same

def quaternion_multiply(q1, q2):

    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2

    return np.array([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2
    ])


def quaternion_conjugate(q):

    return np.array([
        q[0],
        -q[1],
        -q[2],
        -q[3]
    ])


def quaternion_normalize(q):

    return q / np.linalg.norm(q)


def rotate_vector(q, vector):

    vector_quaternion = np.array([
        0.0,
        vector[0],
        vector[1],
        vector[2]
    ])

    q_conjugate = quaternion_conjugate(q)

    rotated = quaternion_multiply(quaternion_multiply(q, vector_quaternion),q_conjugate)

    return rotated[1:]


def inverse_rotate_vector(q, vector):
    return rotate_vector(quaternion_conjugate(q),vector)

def quaternion_to_euler(q):
    w,x,y,z = q

    pitch = np.arctan2(
        2 * (w*x + y*z),
        1 - 2 * (x*x + y*y)
    )

    yaw = np.arctan2(
        2 * (w*y - z*x),
        1 - 2 * (y*y + x*x)
    )

    roll = np.arctan2(
        2 * (w*z + x*y),
        1 - 2 * (z*z + y*y)
    )

    return np.degrees([pitch,yaw,roll])

def quaternion_from_two_vectors(v1, v2):
    v1 = np.array(v1, dtype=float)
    v2 = np.array(v2, dtype=float)

    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)

    dot = np.dot(v1, v2)
    cross = np.cross(v1, v2)

    # Vectors point in the same direction
    if dot > 0.999999:
        return np.array([1.0, 0.0, 0.0, 0.0])

    # Vectors point in opposite directions
    if dot < -0.999999:
        axis = np.cross(v1, np.array([1.0, 0.0, 0.0]))

        if np.linalg.norm(axis) < 0.000001:
            axis = np.cross(v1, np.array([0.0, 1.0, 0.0]))

        axis = axis / np.linalg.norm(axis)

        return np.array([
            0.0,
            axis[0],
            axis[1],
            axis[2]
        ])

    quaternion = np.array([
        1.0 + dot,
        cross[0],
        cross[1],
        cross[2]
    ])

    return quaternion_normalize(quaternion)


position = np.array([0.0, 0.0, 0.0]) #North, East, Up, is world orientated
velocity = np.array([0.0, 0.0, 0.0]) #North, East, Up, is world orientated
orientation = quaternion_from_two_vectors(np.array([0.0,0.0,1.0]),launch_vector) #w, x, y, z. where x is pitch, y is yaw, z is roll, is body orientated
angular_velocity = np.array([0.0, 0.0, 0.0]) #x, y, z, is body orientated
state = [position, velocity, orientation, angular_velocity] #2d array

force = np.array([0.0, 0.0, 0.0]) #North, East, Up, is world orientated
total_moment_body = np.array([0.0, 0.0, 0.0])
momentOfInertia = np.diag([0.0, 0.0, 0.0]) #x, y, z, is body orientated

positionHistory = np.empty((3,0))
velocityHistory = np.empty((3,0))
forceHistory = np.empty((3,0))
orientationHistory = np.empty((3,0))
angularVelocityHistory = np.empty((3,0))
momentHistory = np.empty((3,0))


#Constants
gravity = 9.81 # m/s^2
timeStep = 0.05 # seconds

#calculated values

CT = (CG_wet*(dryMass+fuelMass)-CG_dry*dryMass)/(fuelMass)
CP = 1.149
CN_alpha = 11.97 #pulled from OR

#fileName = input("Enter the CSV file containing motor data: ") # CSV file containing motor data

#motor_csv = import_Motor_Data(fileName)
motor_csv = import_Motor_Data("Motors/Hypertek_L550.csv")
timeCurve = np.array(motor_csv[0])
thrustCurve = np.array(motor_csv[1])

_segment_impulse = 0.5 * (thrustCurve[1:] + thrustCurve[:-1]) * np.diff(timeCurve)
cumulativeImpulse = np.concatenate(([0.0], np.cumsum(_segment_impulse)))
netImpulse = cumulativeImpulse[-1]

#Init Variables
time = 0
leftRod = False

#Init Lists
timeList = []
massList = []
machList = []

def on_rod(position):
    global leftRod
    if leftRod:
        return False
    if np.linalg.norm(position) <= rodLength:
        return True
    leftRod = True
    return False

def thrust_magnitude(time):
    return float(np.interp(time, timeCurve, thrustCurve, right=0.0))

def drag_magnitude(position, velocity):
    A, cd = 3.14 * 0.0508 ** 2, dragCoefficient        
    return (0.5 * 1.1 * cd * A * speed(velocity)**2)

def speed(velocity):
    return np.linalg.norm(velocity)

def weight(time):
    return mass(time) * gravity

def mass(time):
    if time <= timeCurve[0]:
        currentImpulse = thrust_magnitude(time) * time
    elif time >= timeCurve[-1]:
        currentImpulse = netImpulse
    else:
        i = np.searchsorted(timeCurve, time)
        if timeCurve[i] == time:
            currentImpulse = cumulativeImpulse[i]
        else:
            i -= 1
            currentImpulse = cumulativeImpulse[i] + 0.5 * (thrust_magnitude(time) + thrustCurve[i]) * (time - timeCurve[i])

    currentFuelMass = fuelMass * (1 - (currentImpulse / netImpulse))
    return dryMass + currentFuelMass

def centerOfGravity(time):
    return (CG_dry*dryMass+CT*(mass(time)-dryMass))/(mass(time))


def angle_of_attack(velocity):
    if velocity[0] == 0 and velocity[1] == 0:
        return 0.0
    else:
        return np.arctan2(np.linalg.norm(velocity[:2]), velocity[2]) 

def normal_force(velocity):
    dynamicPressure = 0.5*1.1*speed(velocity)**2 #Assuming incompressible air for now
    Area = np.pi * 0.0508 **2
    groundVelocity = np.array([velocity[0], velocity[1], 0.0])
    groundSpeed = np.linalg.norm(groundVelocity)
    if groundSpeed == 0:
        return np.zeros(3)
    direction = -groundVelocity / groundSpeed
    alpha = angle_of_attack(velocity)
    magnitude = (dynamicPressure * Area * CN_alpha * alpha)
    return magnitude * direction

def moment_of_inertia(time):
    # Assuming a simple cylindrical rocket for now
    radius = 0.0508  # m
    length = 1.2  # m
    mass_rocket = mass(time)
    
    I_xx = (1/12) * mass_rocket * (3 * radius**2 + length**2)
    I_yy = I_xx
    I_zz = (1/2) * mass_rocket * radius**2
    
    return ([I_xx, I_yy, I_zz])

def derivatives(time, position, velocity, orientation, angular_velocity):
    wind = np.array([0.0, 0.0, 0.0]) # North, East, Up, is world orientated
    rho = 1.1
    CG = centerOfGravity(time)
    
    # Air-relative velocity in WORLD coordinates
    airFlow_velocity_world = velocity - wind

    # Convert to BODY coordinates
    airFlow_velocity_body = inverse_rotate_vector(orientation, airFlow_velocity_world)

    airFlow_speed = np.linalg.norm(airFlow_velocity_body)

    # Aerodynamics
    drag = drag_magnitude(position, airFlow_velocity_body)
    
    if airFlow_speed != 0:
        drag_body = (-drag * airFlow_velocity_body / airFlow_speed) 
    else:
        drag_body = np.zeros(3) 

    drag_world = rotate_vector(orientation, drag_body)

    aero_force_body = normal_force(airFlow_velocity_body)
    
    # Forces
    thrust_body = np.array([0.0, 0.0, thrust_magnitude(time)])
    thrust_world = rotate_vector(orientation, thrust_body)

    gravity_world = np.array([0.0, 0.0, -mass(time) * gravity])

    netForce_world = thrust_world + gravity_world + drag_world

    # Acceleration
    acceleration = netForce_world / mass(time)

    radius = [0.0, 0.0, CP - CG]  # Vector from center of pressure to center of gravity in body coordinates
    # Moments
    moment = np.cross(radius, aero_force_body)

    # Angular acceleration
    I_matrix = np.diag(moment_of_inertia(time))
    angular_acceleration = np.linalg.solve(I_matrix, moment - np.cross(angular_velocity, I_matrix @ angular_velocity))

    # Quaternion derivative
    omega_quaternion = np.array([
        0.0,
        angular_velocity[0],
        angular_velocity[1],
        angular_velocity[2]
    ])

    ori_dt = 0.5 * quaternion_multiply(orientation,omega_quaternion)

    return (velocity, acceleration, ori_dt, angular_acceleration)

def rk4_step(time, position, velocity, orientation, angular_velocity):

    k1_pos, k1_vel, k1_q, k1_omega = derivatives(time, position, velocity, orientation, angular_velocity)

    k2_pos, k2_vel, k2_q, k2_omega = derivatives(time + timeStep/2, position + k1_pos * timeStep/2, velocity + k1_vel * timeStep/2, orientation + k1_q * timeStep/2, angular_velocity + k1_omega * timeStep/2)

    k3_pos, k3_vel, k3_q, k3_omega = derivatives(time + timeStep/2, position + k2_pos * timeStep/2, velocity + k2_vel * timeStep/2, orientation + k2_q * timeStep/2, angular_velocity + k2_omega * timeStep/2)

    k4_pos, k4_vel, k4_q, k4_omega = derivatives(time + timeStep, position + k3_pos * timeStep, velocity + k3_vel * timeStep, orientation + k3_q * timeStep, angular_velocity + k3_omega * timeStep)


    new_position = position + (timeStep/6) * (k1_pos + 2*k2_pos + 2*k3_pos + k4_pos)
    new_velocity = velocity + (timeStep/6) * (k1_vel + 2*k2_vel + 2*k3_vel + k4_vel)
    new_orientation = orientation + (timeStep/6) * (k1_q + 2*k2_q + 2*k3_q + k4_q)
    new_angular_velocity = angular_velocity + (timeStep/6) * (k1_omega + 2*k2_omega + 2*k3_omega + k4_omega)

    new_orientation = quaternion_normalize(new_orientation)

    if(on_rod(new_position)):
            new_orientation = orientation
            new_angular_velocity = np.zeros(3)
    
    return (new_position,new_velocity,new_orientation,new_angular_velocity)

def run_simulation():
    global time, position, velocity, orientation, angular_velocity, positionHistory, velocityHistory, forceHistory, orientationHistory, angularVelocityHistory, momentHistory


    while -0.01 <= position[2] <= 10000000:  #limits the sim to actual atmosphere, avoid errors.

        euler_orientation = quaternion_to_euler(orientation)

        derivative = derivatives(time, position, velocity, orientation, angular_velocity)
        timeList.append(time)
        massList.append(mass(time))
        positionHistory = np.hstack((positionHistory, np.array(position).reshape(-1, 1)))
        velocityHistory = np.hstack((velocityHistory, np.array(velocity).reshape(-1, 1)))
        forceHistory = np.hstack((forceHistory, np.array(derivative[1] * mass(time)).reshape(-1, 1)))  
        orientationHistory = np.hstack((orientationHistory, np.array(euler_orientation).reshape(-1,1)))
        angularVelocityHistory = np.hstack((angularVelocityHistory, np.array(angular_velocity).reshape(-1, 1)))
        momentHistory = np.hstack((momentHistory, np.array(derivative[3] * np.array(moment_of_inertia(time))).reshape(-1, 1)))

        position, velocity, orientation, angular_velocity = rk4_step(time, position, velocity, orientation, angular_velocity)
        
        time += timeStep


def plot_results():
    fig, ax1 = plt.subplots(figsize=(10, 8))
    ax1.plot(timeList, positionHistory[2], color='blue', label='Altitude (m)')
    ax1.set_ylabel('Altitude (m)', color='blue')
    ax1.tick_params(axis='y', labelcolor='blue')

    ax2 = ax1.twinx()

    ax2.plot(timeList, velocityHistory[2], color='red', label='Velocity (m/s)')
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
        print(f"Time: {timeList[i]:.3f} s, \
              Mass: {massList[i]:.3f} kg, \
              x_pos: {positionHistory[0][i]:.3f} m, \
              y_pos: {positionHistory[1][i]:.3f} m, \
              Altitude: {positionHistory[2][i]:.3f} m, \
              x_vel: {velocityHistory[0][i]:.3f} m/s, \
              y_vel: {velocityHistory[1][i]:.3f} m/s, \
              z_vel: {velocityHistory[2][i]:.3f} m/s, \
              x_force: {forceHistory[0][i]:.3f} N, \
              y_force: {forceHistory[1][i]:.3f} N, \
              z_force: {forceHistory[2][i]:.3f} N, \
              pitch: {orientationHistory[0][i]:.3f} degrees, \
              yaw: {orientationHistory[1][i]:.3f} degrees, \
              roll: {orientationHistory[2][i]:.3f} degrees, \
              pitch_rate: {angularVelocityHistory[0][i]:.3f} rad/s, \
              yaw_rate: {angularVelocityHistory[1][i]:.3f} rad/s, \
              roll_rate: {angularVelocityHistory[2][i]:.3f} rad/s, \
              pitch_moment: {momentHistory[0][i]:.3f} kg*m^2, \
              yaw_moment: {momentHistory[1][i]:.3f} kg*m^2, \
              roll_moment: {momentHistory[2][i]:.3f} kg*m^2")
    rodVelocity = 0
    for i in range(len(positionHistory[0])-1):
        dist1=np.linalg.norm(positionHistory[:, i])
        dist2=np.linalg.norm(positionHistory[:, i+1])
        if dist1 < rodLength <= dist2:
            rodVelocity = np.linalg.norm(velocityHistory[:, i]+velocityHistory[:, i+1])/2
            break
    print(f"Velocity off rod: {rodVelocity:.3f} m/s, Apogee: {max(positionHistory[2]):.3f} m, Max Velocity: {max(velocityHistory[2]):.3f} m/s, Time to Apogee: {timeList[np.argmax(positionHistory[2])]:.3f} s, Flight Time: {timeList[-1]:.3f} s, Ground hit velocity: {velocityHistory[2][-1]:.3f} m/s, Range from launch: {np.sqrt(positionHistory[0][-1]**2 + positionHistory[1][-1]**2):.3f} m")
    sys.stdout.flush()
    plot_results()
