#sim.py - V2.3
import matplotlib.pyplot as plt 
import numpy as np
import csv
import sys
import random
import json

parachuteDragCoefficient = 0  # not stored numerically ("auto" in file)
parachuteArea = 0
finSets = 1
transitions = 0
MainDeploymentAltitude = 0
MainArea = 0
MainDragCoefficient = 0
epsilonAngle = 0 #basically transition angle for nosecone, 0 since smooth, but a cone could be like 15
surfaceRoughness = 0    # meters — needs a materials-roughness table/guess
launchGuideRoughness = 0    # meters — same as above
dryMass = 0  # kg — mass properties computed during simulation
fuelMass = 0   # kg — motor propellant mass
rodLength   = 1                # m — launch condition, not part of <rocket>
launchAngle = 0
launchDirection = 0 #heading
launch_vector = np.array([0.0, 0.0, 0.0])
CG_dry = 0  # m from nose tip — computed mass property
CG_wet = 0  # m from nose tip — same
leadingEdgeToNosecone = 0
bodyLength=0
motorDiameter = 0
launchGuideLength, launchGuideOuterDiameter, launchGuideInnerDiameter = 0,0,0
rodVelocity = None
cpBiasZ = 0.0          # per-run draw ~N(0,1); 0 = nominal
cpSigmaBase = 0.10     # cal, 1-sigma away from Mach 1
cpSigmaPeak = 0.30     # cal, extra 1-sigma at Mach 1
cpWidth = 0.22         # Mach half-width of the bump
cpJitter = 0.0         # cal, time-varying wander (0 = off)
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
orientation = np.array([0.0, 0.0, 0.0, 0.0]) #w, x, y, z. where x is pitch, y is yaw, z is roll, is body orientated
angular_velocity = np.array([0.0, 0.0, 0.0]) #x, y, z, is body orientated

positionHistory = np.empty((3,0))
velocityHistory = np.empty((3,0))
forceHistory = np.empty((3,0))
orientationHistory = np.empty((3,0))
angularVelocityHistory = np.empty((3,0))
momentHistory = np.empty((3,0))
numberOfFins = 4

#Constants
gravity = 9.806 # m/s^2
timeStep = 0.05 # seconds
airMolarMass = 0.0289644 #kg/mol
universeGasConstant = 8.31447 #J/(mol*K)
SpecificGasConstant = 287.05 #J/(kg*K)
temperatureLapseRate = 0.0065 #K/m
AdiabaticIndex = 1.4 #unitless & constant for our purpose, value for air
initialTemperature = 0
initialAltitude = 0
initialPressure = 0
avgWindSpeed = 0
turbulence = 0 
windHeading = 0
x_n1 = 0
x_n2 = 0
finSpan, rootChord, tipChord, sweepLength, noseLength, bodytubeLength, finThickness = 0,0,0,0,0,0,0
bodyDiameter = 0
#calculated values
Y1 = 0 #supposed to be lambda or the distance on Fin at the mid-chord lines for set 1 (m)
A_fin = 0
finenessRatio = 0 #unitless, ratio of the length of the rocket to the diameter of the rocket
meanAerodynamicChordLengthOfFin = 0
A_wet_nose = 0
noseVolume = 0
A_wet_body = 0
A_wet_fins = 0
A_ref = 0
CP = 0

cpNoise = []

cpNoiseTimes = [i*0.1 for i in range(10000)]

def initcenterOfPressure():
    CN_n = 2
    CN_f1 = (1+((bodyDiameter/2))/(finSpan+(bodyDiameter/2))) * ((4*numberOfFins*(finSpan/bodyDiameter)**2)/(1+np.sqrt(1+((2*Y1)/(rootChord+tipChord))**2)))
    Pf1 = CN_f1 * (leadingEdgeToNosecone+(sweepLength*(rootChord+2*tipChord)/(3*(rootChord+tipChord)))+(1/6)*(rootChord+tipChord-(rootChord*tipChord)/(rootChord+tipChord)))
    Pn = CN_n * (0.466*noseLength)
    

    netNormalForce = CN_n + CN_f1
    netMoment = Pn + Pf1
    return netMoment/netNormalForce

#motor_csv = import_Motor_Data(fileName)
timeCurve = np.empty(0)
thrustCurve = np.empty(0)

_segment_impulse = 0
cumulativeImpulse = 0
netImpulse = 0
#Init Variables
time = 0
leftRod = False

#Init Lists
timeList = []
massList = []
machList = []
AOAList = []
CDList = []

timeWindSpeeds = []
windSpeeds = []

def temperature(position):
    altitude = position[2]
    if altitude+initialAltitude < 11000:
        return initialTemperature - temperatureLapseRate * (altitude)
    else:
        return initialTemperature - temperatureLapseRate * (11000-initialAltitude)

def pressure(position):
    return initialPressure * ((temperature(position)+273.15)/(initialTemperature+273.15))**(gravity*airMolarMass/(universeGasConstant*temperatureLapseRate))

def airDensity(position):
    return pressure(position)/(SpecificGasConstant*(temperature(position)+273.15))

def speedOfSound(position):
    return np.sqrt(AdiabaticIndex*SpecificGasConstant*(temperature(position)+273.15))

def machNumber(airVelocity, position):
    return speed(airVelocity)/speedOfSound(position)

def speed(velocity):
    return np.linalg.norm(velocity)

def angle_of_attack(velocity):
    if velocity[0] == 0 and velocity[1] == 0:
        return 0.0
    else:
        return np.arctan2(np.linalg.norm(velocity[:2]), velocity[2]) 
    
def dynamicViscosity(position):
    return 1.458*10**(-6)*((temperature(position)+273.15)**(3/2))/(temperature(position)+273.15+110.4)

def kinematicViscosity(position):
    return dynamicViscosity(position)/airDensity(position)

def reynoldsNumber(airVelocity, position, characteristicLength):
    return speed(airVelocity) * characteristicLength / kinematicViscosity(position)

def zeroAngleDragCoefficient(time, airVelocity, position):
    C_d = C_D_friction(airVelocity, position) + nosePressureDrag(airVelocity, position) + finPressureDrag(time, airVelocity, position) + baseDragCoefficient(time, airVelocity, position) + parasiticDrag(airVelocity, position)
    return C_d

def C_D_friction(airVelocity, position):
    R = reynoldsNumber(airVelocity, position, (bodytubeLength+noseLength))
    M = machNumber(airVelocity, position)
    R_crit = 51*(surfaceRoughness/(bodytubeLength+noseLength))**-1.039
    if R < 10**4:
        C_f = 1.48*10**-2
    elif R < R_crit:
        C_f = 1/(1.5*np.log(R)-5.6)**2
    else:
        C_f = 0.032*(surfaceRoughness/(bodytubeLength+noseLength))**0.2

    C_Mach = (1-0.1*M**2)
    C_f_rough = 0.032*(surfaceRoughness/(bodytubeLength+noseLength))**0.2 * C_Mach
    C_f_component = max(C_f, C_f_rough)
    K_body = 1 + 1/(2*finenessRatio)
    K_fin = 1 + 2*finThickness/meanAerodynamicChordLengthOfFin
    C_D_body = C_f_component * K_body * (A_wet_body/A_ref)
    C_D_fins = C_f_component * K_fin * (A_wet_fins/A_ref)
    return C_D_body + C_D_fins

def nosePressureDrag(airVelocity, position):
    jointAngle_phi=0.0
    kappa=1.0
    M = machNumber(airVelocity, position)
    gamma = 1.4
    eps = np.arctan(bodyDiameter/(2*noseLength))          # (B.3)
    s = np.sin(eps)
    K = 0.72*(kappa-0.5)**2 + 0.82                        # (B.8)
    C0 = 0.8*np.sin(jointAngle_phi)**2                    # (3.86)
    C1 = K*s                                              # (B.6)
    dC1 = K*4/(gamma+1)*(1-0.5*s)                         # (B.5)
    def sup(m):  return K*(2.1*s**2 + 0.5*s/np.sqrt(m**2-1))   # (B.4)
    def dsup(m): return -K*0.5*s*m/(m**2-1)**1.5
    if M <= 1:                                            # (3.87)
        a = C1 - C0
        return a*M**(dC1/a) + C0
    if M < 1.3:                                           # cubic Hermite, 1 -> 1.3
        h = 0.3; t = (M-1)/h
        return ((2*t**3-3*t**2+1)*C1 + (t**3-2*t**2+t)*dC1*h
                + (-2*t**3+3*t**2)*sup(1.3) + (t**3-t**2)*dsup(1.3)*h)
    return sup(M)  

def finPressureDrag(time, airVelocity, position):
    M = machNumber(airVelocity, position)
    if M < 0.9: #rounded edge
        perpindicularLeadingEdgeDrag = (1-M**2)**(-0.417) -1
    elif M < 1:
        perpindicularLeadingEdgeDrag = 1 - 1.785*(M-0.9)
    else:
        perpindicularLeadingEdgeDrag = 1.214 - (0.502/M**2) + (0.1095/M**4)
    leadingEdgeDrag = perpindicularLeadingEdgeDrag * np.cos(leadingEdgeAngle)**2 #angled fin, sweep angle
    trailingEdgeDrag = (1/2) * baseDragCoefficient(time, airVelocity, position) #rounded edge
    C_D_fin = leadingEdgeDrag + trailingEdgeDrag
    A_FIN = numberOfFins * finThickness * finSpan
    return (A_FIN/A_ref) * C_D_fin

def baseDragCoefficient(time, airVelocity, position):
    if time < timeCurve[-1]:
        A_motor = (np.pi/4) * motorDiameter ** 2
    else:
        A_motor = 0
    M = machNumber(airVelocity, position)
    if M < 1:
        CD_base = 0.12 + 0.13*M**2
    else:
        CD_base = 0.25/M
    return CD_base * (A_ref-A_motor)/A_ref

def parasiticDrag(airVelocity, position):
    M = machNumber(airVelocity, position)
    if M < 1:
        qRatio = 1 + (M**2/4) + (M**4/40)
        C_D_base = 0.12 + 0.13*M**2
    else:
        qRatio = 1.84 - (0.76/M**2) + (0.166/M**4) + (0.035/M**6)
        C_D_base = 0.25/M
    C_D_stag = 0.85 * qRatio
    C_D_parasitic = max(1.3-0.3*(launchGuideLength/launchGuideOuterDiameter), 1) * C_D_stag
    outer_Area = np.pi * (launchGuideOuterDiameter/2)**2
    inner_Area = np.pi * (launchGuideInnerDiameter/2)**2
    A_parasitic = outer_Area - inner_Area
    if np.linalg.norm(airVelocity) == 0:
        return 0
    else:
        f = 0.25/((np.log10((launchGuideRoughness/(3.7*launchGuideInnerDiameter)) + (5.74/(reynoldsNumber(airVelocity, position, launchGuideInnerDiameter)**0.9))))**2)
    C_D_tube = f * (launchGuideLength/launchGuideInnerDiameter)
    return (C_D_tube * inner_Area + 0.7*(C_D_parasitic + C_D_base)*A_parasitic)/A_ref

def axialDragCoefficient(time, airVelocity, position, AOA):
    if abs(np.degrees(AOA)) <= 17:
        scalefunction = (-3/24565)*abs(np.degrees(AOA))**3 + (9/2890)*(np.degrees(AOA))**2 + 1
    else:
        scalefunction = (13/1945085)*(np.degrees(AOA)-90)**3 + (2847/3890170)*(np.degrees(AOA)-90)**2
    return zeroAngleDragCoefficient(time, airVelocity, position) * scalefunction

def axialDragForce(time, airVelocity, position, AOA):
    rho = airDensity(position)
    V = speed(airVelocity)
    if airVelocity[2] >= 0:
        A, C_A = A_ref, axialDragCoefficient(time, airVelocity, position, AOA)
    elif position[2] > MainDeploymentAltitude:
        A, C_A = parachuteArea, parachuteDragCoefficient
    else:
        A, C_A = MainArea, MainDragCoefficient
    return (1/2)*rho*V**2 * A * C_A

def normalForceFinCoefficientDerivative(airVelocity, position):
    M = machNumber(airVelocity, position)
    midchordSweepAngle = np.arctan((sweepLength + tipChord/2 - rootChord/2) / finSpan) if finSpan else 0.0
    beta = np.sqrt(abs(1-M**2))
    C_Nalpha0 = 2*np.pi/beta
    AR = 2*(finSpan**2)/A_fin
    F_D = AR/((1/(2*np.pi))*C_Nalpha0*np.cos(midchordSweepAngle))
    C_N_1fin = (C_Nalpha0 * F_D * (A_fin/A_ref) * np.cos(midchordSweepAngle))/(2+F_D*np.sqrt(1+4/(F_D**2)))    #(2*np.pi*((finSpan**2)/A_ref))/(2+((beta*finSpan**2)/(A_fin*np.cos(midchordSweepAngle))))
    C_N_fins = (numberOfFins/2) * C_N_1fin
    K_TB = 1 + ((bodyDiameter/2)/(finSpan+(bodyDiameter/2)))
    return K_TB * C_N_fins


def normalForceCoefficientDerivative(airVelocity, position,AOA):
    C_N_finsWithInterference = normalForceFinCoefficientDerivative(airVelocity, position)
    C_N_nose = 2*np.cos(AOA)
    return C_N_nose + C_N_finsWithInterference

def pitchMomentCoefficientDerivative(airVelocity, position, AOA):
    C_m_nose = (2*np.sin(AOA))/(A_ref*bodyDiameter)*(noseLength*A_ref-noseVolume)
    M = machNumber(airVelocity, position)
    AR = 2*finSpan**2/A_fin
    beta = np.sqrt(abs(1-M**2))

    A = np.array([
    [1,0.5,0.5**2,0.5**3,0.5**4,0.5**5],
    [0,1,2*0.5,3*0.5**2,4*0.5**3,5*0.5**4],
    [1,2,2**2,2**3,2**4,2**5],
    [0,1,2*2,3*2**2,4*2**3,5*2**4],
    [0,0,2,6*2,12*2**2, 20*2**3],
    [0,0,0,6,24*2,60*2**2]
    ])

    b = np.array([
        0.25, #p(0.5)
        0,  #p'(0.5)
        meanAerodynamicChordLengthOfFin*((AR*np.sqrt(abs(2**2-1)) - 0.67)/(2*AR*np.sqrt(abs(2**2-1))-1)),   #p(2) = f(2)
        ((0.68*AR)/(np.sqrt(3)*(2*AR*np.sqrt(3)-1)**2)),  #p'(2) = f'(2)
        0,  #p''(2)
        0   #p'''(2)
    ])

    coefficients = np.linalg.solve(A,b)

    if M <= 0.5:
        x_f = (sweepLength/3)*((rootChord+2*tipChord)/(rootChord+tipChord))+(1/6)*((rootChord**2+tipChord**2+rootChord*tipChord)/(rootChord+tipChord))
    elif M <=2:
        x_f = coefficients[0] + coefficients[1]*M + coefficients[2]*M**2 + coefficients[3]*M**3 + coefficients[4]*M**4 + coefficients[5]*M**5
    else:
        x_f = meanAerodynamicChordLengthOfFin*((AR*beta - 0.67)/(2*AR*beta-1))
    x_fin = leadingEdgeToNosecone + x_f
    C_m_fins = normalForceFinCoefficientDerivative(airVelocity, position) * x_fin/bodyDiameter
    return C_m_nose+C_m_fins

def normalForce(airVelocity, position, AOA):
    rho = airDensity(position)
    V = speed(airVelocity)
    C_N_alpha = normalForceCoefficientDerivative(airVelocity, position, AOA)
    return (1/2)*rho*V**2 * A_ref * C_N_alpha * AOA

def centerOfPressure(airVelocity, position, AOA):
    if AOA != 0:
        return (pitchMomentCoefficientDerivative(airVelocity, position, AOA)/normalForceCoefficientDerivative(airVelocity, position, AOA))*bodyDiameter
    else:
        return CP

def moment_of_inertia(time):
    # Assuming a simple cylindrical rocket for now
    Mass = mass(time)
    CG = centerOfGravity(time)
    
    I_xx = Mass * CG * ((bodytubeLength+noseLength)-CG)
    I_yy = I_xx
    I_zz = Mass * (bodyDiameter/2)**2
    
    return ([I_xx, I_yy, I_zz])

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

# --- CP uncertainty (calibers) ---

def cpSigma(M):
    return cpSigmaBase + cpSigmaPeak*np.exp(-((M-1)/cpWidth)**2)

def cpOffset(time, M):
    """CP shift in metres (+ = aft = more stable)."""
    bump = np.exp(-((M-1)/cpWidth)**2)
    wander = float(np.interp(time, cpNoiseTimes, cpNoise)) if cpJitter else 0.0
    return (cpBiasZ*cpSigma(M) + cpJitter*bump*wander) * bodyDiameter



def derivatives(time, position, velocity, orientation, angular_velocity):
    wind = float(np.interp(time, timeWindSpeeds, windSpeeds, right=0.0))
    windVector = wind*np.array([np.sin(windHeading), np.cos(windHeading), 0])
    airFlow_velocity_world = velocity - windVector
    airFlow_velocity_body = inverse_rotate_vector(orientation, airFlow_velocity_world)
    AOA = angle_of_attack(airFlow_velocity_body)
    CP = centerOfPressure(airFlow_velocity_world, position, AOA) + cpOffset(time, machNumber(airFlow_velocity_world, position))
    CG = centerOfGravity(time)

    airFlow_speed = np.linalg.norm(airFlow_velocity_body)


    perpindicular_airFlow = airFlow_velocity_world - airFlow_speed * np.cos(AOA) * rotate_vector(orientation, [0,0,1])
    if AOA != 0:
        unit_perpindicular_airFlow = (1/np.linalg.norm(perpindicular_airFlow)) * perpindicular_airFlow
    else:
        unit_perpindicular_airFlow = np.zeros(3)

    normal_Force_body = inverse_rotate_vector(orientation, -normalForce(airFlow_velocity_world, position, AOA)*unit_perpindicular_airFlow)

    axialDrag_Force = axialDragForce(time, airFlow_velocity_world, position, AOA)
    if airFlow_velocity_world[2] < 0 and airFlow_speed > 0:
        aero_force_world = -axialDrag_Force * airFlow_velocity_world / airFlow_speed
        aero_force_body = inverse_rotate_vector(orientation, aero_force_world)
        moment_scale = 0.0
    else:
        aero_force_body = normal_Force_body + (np.array([0,0,axialDrag_Force]) * -np.sign(airFlow_velocity_body[2]))
        aero_force_world = rotate_vector(orientation, aero_force_body)
        moment_scale = 1.0
    thrust_body = np.array([0.0, 0.0, thrust_magnitude(time)])
    thrust_world = rotate_vector(orientation, thrust_body)

    gravity_world = np.array([0.0, 0.0, -mass(time) * gravity])
    if on_rod(position):
        netForce_world = thrust_world
    else:
        netForce_world = thrust_world + gravity_world + aero_force_world

    # Acceleration
    acceleration = netForce_world / mass(time)

    # Moments
    total_moment_body = moment_scale * np.cross([0,0,CG-CP], aero_force_body)

    # Angular acceleration
    I_matrix = np.diag(moment_of_inertia(time))
    angular_acceleration = np.linalg.solve(I_matrix, total_moment_body - np.cross(angular_velocity, I_matrix @ angular_velocity))

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

def run_simulation(INITTEMP, INITALT, INITPRE, AVGSPEED, TURB, HEAD, PARACD, PARAAREA, DEPALT , MAINAREA, MAINCD, MAINROUGH, GUIDEROUGH, DRYM, FUELM, LAUNCHANG, LAUNCHDIR, CGDRY, CGWET, BDIA, THRUSTC, TIMEC, FSPAN, ROOT, TIP, SWEEP, NOSELENGTH, BODYLENGTH, THICK, MOTORDIA, GUIDEL, GUIDEOD, GUIDEID, FINNUM):
    global time, position, velocity, orientation, angular_velocity, positionHistory, velocityHistory, forceHistory, orientationHistory, angularVelocityHistory, momentHistory, machList, initialTemperature, initialAltitude, initialPressure, avgWindSpeed, turbulence, windHeading, parachuteDragCoefficient, parachuteArea, MainDeploymentAltitude, MainArea, MainDragCoefficient, surfaceRoughness, launchGuideRoughness, dryMass, fuelMass, launchAngle, launchDirection, CG_dry, CG_wet, CT, bodyDiameter, thrustCurve, timeCurve, _segment_impulse, cumulativeImpulse, netImpulse,  finSpan, rootChord, tipChord, sweepLength, noseLength, bodytubeLength, leadingEdgeToNosecone, finThickness, Y1, A_fin, finenessRatio, meanAerodynamicChordLengthOfFin, A_wet_nose, noseVolume, A_wet_body, A_wet_fins, A_ref, leadingEdgeAngle, motorDiameter, launchGuideLength, launchGuideOuterDiameter, launchGuideInnerDiameter, CP, rodVelocity, leftRod, timeList, massList, machList, AOAList, CDList, timeWindSpeeds, windSpeeds, x_n1, x_n2, numberOfFins, cpNoise

    initialTemperature, initialAltitude, initialPressure, avgWindSpeed, turbulence, windHeading, parachuteDragCoefficient, parachuteArea, MainDeploymentAltitude, MainArea, MainDragCoefficient, surfaceRoughness, launchGuideRoughness, dryMass, fuelMass, launchAngle, launchDirection, CG_dry, CG_wet, bodyDiameter, thrustCurve, timeCurve, finSpan, rootChord, tipChord, sweepLength, noseLength, bodytubeLength, finThickness, motorDiameter, launchGuideLength, launchGuideOuterDiameter, launchGuideInnerDiameter, numberOfFins = INITTEMP, INITALT, INITPRE, AVGSPEED, TURB, HEAD, PARACD, PARAAREA, DEPALT , MAINAREA, MAINCD, MAINROUGH, GUIDEROUGH, DRYM, FUELM, LAUNCHANG, LAUNCHDIR, CGDRY, CGWET, BDIA, THRUSTC, TIMEC, FSPAN, ROOT, TIP, SWEEP, NOSELENGTH, BODYLENGTH, THICK, MOTORDIA, GUIDEL, GUIDEOD, GUIDEID, FINNUM
    CT = (CG_wet*(dryMass+fuelMass)-CG_dry*dryMass)/(fuelMass)
    _segment_impulse = 0.5 * (thrustCurve[1:] + thrustCurve[:-1]) * np.diff(timeCurve)
    cumulativeImpulse = np.concatenate(([0.0], np.cumsum(_segment_impulse)))      
    netImpulse = cumulativeImpulse[-1]
    leadingEdgeToNosecone = (bodytubeLength+noseLength) - rootChord

    launch_vector = np.array([np.sin(launchAngle)*np.cos(launchDirection), np.sin(launchAngle) * np.sin(launchDirection), np.cos(launchAngle)])

    Y1 = np.sqrt(finSpan**2+(sweepLength+tipChord/2-rootChord/2)**2) #supposed to be lambda or the distance on Fin at the mid-chord lines for set 1 (m)
    A_fin = (finSpan/2)*(tipChord+rootChord)
    finenessRatio = (bodytubeLength+noseLength)/bodyDiameter #unitless, ratio of the length of the rocket to the diameter of the rocket
    meanAerodynamicChordLengthOfFin = (2/3) * (rootChord + tipChord - (rootChord*tipChord)/(rootChord + tipChord)) #meters, this is a typical value for a fin
    A_wet_nose = (np.pi/(2*(bodyDiameter/2)**2))*((bodyDiameter/2)*noseLength*((bodyDiameter/2)**2-noseLength**2)+(((bodyDiameter/2)**2+noseLength**2)**2)*np.arctan((bodyDiameter/2)/noseLength)) #meters^2, this is the wetted area of an ogive nosecone of the rocket
    noseVolume = (np.pi/(24*(bodyDiameter/2)**3))*(6*noseLength*(bodyDiameter/2)**5+6*(bodyDiameter/2)*noseLength**5+4*((bodyDiameter/2)**3)*(noseLength**3)+6*((bodyDiameter/2)**2-noseLength**2)*((bodyDiameter/2)**2+noseLength**2)**2*np.arctan((bodyDiameter/2)/noseLength))
    A_wet_body = np.pi * bodyDiameter * bodytubeLength + A_wet_nose #meters^2, this is the wetted area of the body of the rocket
    A_wet_fins = numberOfFins*(finSpan*(rootChord+tipChord)+finThickness*(rootChord+tipChord+np.sqrt(sweepLength**2 + finSpan**2)+np.sqrt((sweepLength+tipChord-rootChord)**2+finSpan**2))) #meters^2, this is the wetted area of the fins of the rocket
    A_ref = np.pi * ((bodyDiameter/2)**2) #meters^2, this is the reference area of the rocket, which is the cross-sectional area of the rocket body

    CP = initcenterOfPressure() # m from the nose tip

    leadingEdgeAngle = np.arctan(sweepLength / finSpan)

    time = 0
    position = np.array([0.0, 0.0, 0.0]) #North, East, Up, is world orientated
    velocity = np.array([0.0, 0.0, 0.0]) #North, East, Up, is world orientated
    orientation = quaternion_from_two_vectors(np.array([0.0,0.0,1.0]),launch_vector) #w, x, y, z. where x is pitch, y is yaw, z is roll, is body orientated



    positionHistory = np.empty((3,0))
    velocityHistory = np.empty((3,0))
    forceHistory = np.empty((3,0))
    orientationHistory = np.empty((3,0))
    angularVelocityHistory = np.empty((3,0))
    momentHistory = np.empty((3,0))
    rodVelocity = None
    leftRod = False

        
    timeList = []
    massList = []
    machList = []
    AOAList = []
    CDList = []

    timeWindSpeeds = []
    windSpeeds = []
    x_n1 = 0
    x_n2 = 0

   
    c1, c2 = 0.0, 0.0
    for i in range(10000):              # same AR(2) filter as the wind
        c = random.gauss(0,1) + (5/6)*c1 - (5/24)*c2
        c2, c1 = c1, c
        cpNoise.append(c)

    for i in range(0,10000):
        w_n = random.gauss(0,1)
        x_n = w_n + (5/6)*(x_n1)-(5/24)*(x_n2)
        x_n2 = x_n1
        x_n1 = x_n
        stdev = avgWindSpeed*turbulence
        windSpeeds.append(avgWindSpeed + stdev*x_n)
        timeWindSpeeds.append(i*0.1)



    while -0.01 <= position[2] <= 10000000:  #limits the sim to actual atmosphere, avoid errors.

        euler_orientation = quaternion_to_euler(orientation)
        wind = float(np.interp(time, timeWindSpeeds, windSpeeds, right=0.0))
        windVector = wind* np.array([np.sin(windHeading),np.cos(windHeading), 0])
        derivative = derivatives(time, position, velocity, orientation, angular_velocity)
        timeList.append(time)
        massList.append(mass(time))
        positionHistory = np.hstack((positionHistory, np.array(position).reshape(-1, 1)))
        velocityHistory = np.hstack((velocityHistory, np.array(velocity).reshape(-1, 1)))

        forceHistory = np.hstack((forceHistory, np.array(derivative[1] * mass(time)).reshape(-1, 1)))  
        orientationHistory = np.hstack((orientationHistory, np.array(euler_orientation).reshape(-1,1)))
        angularVelocityHistory = np.hstack((angularVelocityHistory, np.array(angular_velocity).reshape(-1, 1)))
        momentHistory = np.hstack((momentHistory, np.array(derivative[3] * np.array(moment_of_inertia(time))).reshape(-1, 1)))
        machList.append(machNumber(velocity - windVector, position))
        AOAList.append(angle_of_attack(velocity))
        if velocity[2] >= 0:
            #CDList.append(parasiticDrag(velocity,position))
            #CDList.append(baseDragCoefficient(time, velocity, position))
            #CDList.append(finPressureDrag(time, velocity, position))
            #CDList.append(nosePressureDrag(velocity, position))
            #CDList.append(C_D_friction(velocity, position))
            CDList.append(axialDragCoefficient(time, velocity - windVector, position, angle_of_attack(velocity)))
        else:
            CDList.append(parachuteDragCoefficient)


        position, velocity, orientation, angular_velocity = rk4_step(time, position, velocity, orientation, angular_velocity)

        time += timeStep
    for i in range(len(positionHistory[0])-1):
        dist1=np.linalg.norm(positionHistory[:, i])
        dist2=np.linalg.norm(positionHistory[:, i+1])
        if dist1 < rodLength <= dist2:
            rodVelocity = np.linalg.norm((velocityHistory[:, i]+velocityHistory[:, i+1])/2)
            break
    if rodVelocity is None:
        rodVelocity = 0
    machTime = 0
    for i in range(len(timeList)-1):
        if machList[i] >= 1:
            machTimei = i
            break
    if machTimei:
        for i in range(machTimei+1,len(timeList)-1):
            if machList[i] <= 1:
                machTimef = i
                break
    machTime = timeList[machTimef] - timeList[machTimei]
    return (rodVelocity, max(positionHistory[2]), max(velocityHistory[2]), timeList[np.argmax(positionHistory[2])], timeList[-1], velocityHistory[2][-1], np.sqrt(positionHistory[0][-1]**2 + positionHistory[1][-1]**2), max(machList), machTime, timeList, positionHistory, velocityHistory)

    #return positionHistory, velocityHistory, machList, massList, forceHistory, orientationHistory, angularVelocityHistory, momentHistory, AOAList, CDList
#(rodVelocity, max(positionHistory[2]), max(velocityHistory[2]), timeList[np.argmax(positionHistory[2])], timeList[-1], velocityHistory[2][-1], np.sqrt(positionHistory[0][-1]**2 + positionHistory[1][-1]**2), max(machList))

