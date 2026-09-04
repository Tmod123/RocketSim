#sim.py - V2.3

import matplotlib.pyplot as plt 
import numpy as np
import csv
import sys

def import_Motor_Data(filename):
    times = [0.000]
    thrusts = [0.000]
    with open(filename, newline='') as f:
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

#Constants
gravity = 9.81 # m/s^2
timeStep = 0.05 # seconds
initialTemperature = 20 #degrees Celsius
initialAltitude = 220 #meters above sea level
initialPressure = 101550 #Pa
airMolarMass = 0.0289644 #kg/mol
universeGasConstant = 8.31447 #J/(mol*K)
SpecificGasConstant = 287.05 #J/(kg*K)
temperatureLapseRate = 0.0065 #K/m
hellMannCoefficient = 0.34 #constant for Neutral air above human-inhibited area
AdiabaticIndex = 1.4 #unitless & constant for our purpose, value for air
initWindSpeed = 0 #m/s
windHeading = 0 #radians, 0 is wind to the north, pi/2 is wind to the east, pi is wind to the south, 3pi/2 is wind to the west
surfaceRoughness = 60*10**-6 #meters, first number is micrometer, you have to look in a table or guesstimate for this value. 
launchGuideRoughness = 60*10**-9
numberOfFins = 3
rocketLength = 0.425 #meters
finThickness = 0.002 #meters
rootChord = 0.0508 #meters
tipChord = 0.0508 #meters
sweepLength = 0.0254 #meters
leadingEdgeAngle = np.radians(40.3)
finSpan = 0.03 #meters
Y1 = np.sqrt(finSpan**2+(sweepLength+tipChord/2-rootChord/2)**2) #supposed to be lambda or the distance on Fin at the mid-chord lines for set 1 (m)
leadingEdgeToNosecone = 0.349 #meters 
bodyDiameter = 0.025 #meters
bodyRadius = bodyDiameter/2 #meters
motorDiameter = 0.018 # meters also just the inner body tube innerdiameter
noseLength = 0.1 #meters
A_fin = (finSpan/2)*(tipChord+rootChord)
finenessRatio = rocketLength/bodyDiameter #unitless, ratio of the length of the rocket to the diameter of the rocket
meanAerodynamicChordLengthOfFin = (2/3) * (rootChord + tipChord - (rootChord*tipChord)/(rootChord + tipChord)) #meters, this is a typical value for a fin
A_wet_nose = (np.pi/(2*bodyRadius**2))*(bodyRadius*noseLength*(bodyRadius**2-noseLength**2)+((bodyRadius**2+noseLength**2)**2)*np.atan(bodyRadius/noseLength)) #meters^2, this is the wetted area of an ogive nosecone of the rocket
noseVolume = (np.pi/(24*bodyRadius**3))*(6*noseLength*bodyRadius**5+6*bodyRadius*noseLength**5+4*(bodyRadius**3)*(noseLength**3)+6*(bodyRadius**2-noseLength**2)*(bodyRadius**2+noseLength**2)**2*np.arctan(bodyRadius/noseLength))
A_wet_body = np.pi * bodyDiameter * (rocketLength-noseLength) + A_wet_nose #meters^2, this is the wetted area of the body of the rocket
A_wet_fins = numberOfFins*(finSpan*(rootChord+tipChord)+finThickness*(rootChord+tipChord+np.sqrt(sweepLength**2 + finSpan**2)+np.sqrt((sweepLength+tipChord-rootChord)**2+finSpan**2))) #meters^2, this is the wetted area of the fins of the rocket
A_ref = np.pi * (bodyRadius**2) #meters^2, this is the reference area of the rocket, which is the cross-sectional area of the rocket body
jointAngle = 0 #radians, this is the angle of the joint between the body and the fin, which is typically 0 for a rocket with fins that are perpendicular to the body
launchGuideLength = 0.035 #meters
launchGuideOuterDiameter = 0.007 # meters
launchGuideInnerDiameter = 0.005 #meters

finSets = 1
transitions = 0
noseType = "ogive" #cone, ogive, paraboloid, ellipsoid

def initcenterOfPressure():
    CN_n = 2
    CN_f1 = (1+(bodyRadius)/(finSpan+bodyRadius)) * ((4*numberOfFins*(finSpan/bodyDiameter)**2)/(1+np.sqrt(1+((2*Y1)/(rootChord+tipChord))**2)))
    Pf1 = CN_f1 * (leadingEdgeToNosecone+(sweepLength*(rootChord+2*tipChord)/(3*(rootChord+tipChord)))+(1/6)*(rootChord+tipChord-(rootChord*tipChord)/(rootChord+tipChord)))

    if noseType == "cone":
        Pn = CN_n * (0.6667*noseLength)
    elif noseType == "ogive":
        Pn = CN_n * (0.466*noseLength)
    elif noseType == "paraboloid":
        Pn = CN_n * (0.5*noseLength)
    elif noseType == "ellipsoid":
        Pn = CN_n * (0.3333*noseLength)
    else:
        Pn = 0

    netNormalForce = CN_n + CN_f1
    netMoment = Pn + Pf1
    return netMoment/netNormalForce

#everyting else
dryMass = 0.0605 # kg
fuelMass = 0.011 # kg
radius = 0.0127 # m
dragCoefficient = 0.634 # dimensionless
parachuteArea = 0.0707 # m^2
parachuteDragCoefficient = 0.80 # dimensionless
rodLength = 1.0 # m
launchAngle = np.radians(5) #radians of x degrees
windSpeed = 0 # m/s
rocketLength = 0.425 # m
CG_dry = 0.24 # m from the nose tip
CG_wet = 0.26
CP = initcenterOfPressure() # m from the nose tip
CT = (CG_wet*(dryMass+fuelMass)-CG_dry*dryMass)/(fuelMass)

CN_alpha = 11.97 #pulled from OR

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
theta = launchAngle
omega = 0
leftRod = False

#Init Lists
timeList = []
massList = []
x_posList = []
y_posList = []
x_velocityList = []
y_velocityList = []
thetaList = []
omegaList = []
alphaList = []
torqueList = []
inertiaList = []
normalForceList = []
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

def temperature(altitude):
    if altitude+initialAltitude < 11000:
        return initialTemperature - temperatureLapseRate * (altitude)
    else:
        return initialTemperature - temperatureLapseRate * (11000-initialAltitude)

def pressure(altitude):
    return initialPressure * ((temperature(altitude)+273.15)/(initialTemperature+273.15))**(gravity*airMolarMass/(universeGasConstant*temperatureLapseRate))

def dynamicPressure(windSpeed, rocketVelocityX, rocketVelocityY, altitude):
    return 0.5 * airDensity(altitude) * airSpeed(windSpeed, rocketVelocityX, rocketVelocityY)**2

def airDensity(altitude):
    return pressure(altitude)/(SpecificGasConstant*(temperature(altitude)+273.15))

def speedOfSound(altitude):
    return np.sqrt(AdiabaticIndex*SpecificGasConstant*(temperature(altitude)+273.15))

def machNumber(windSpeed, rocketVelocityX, rocketVelocityY, altitude):
    return airSpeed(windSpeed, rocketVelocityX, rocketVelocityY)/speedOfSound(altitude)

def airSpeed(windspeed, rocketVelocityX, rocketVelocityY):
    return np.sqrt((rocketVelocityX-windspeed*np.sin(windHeading))**2 + (rocketVelocityY**2) + (windspeed*np.cos(windHeading))**2)

def rocketSpeed(rocketVelocityX, rocketVelocityY):
    return np.sqrt(rocketVelocityX**2 + rocketVelocityY**2)

def angle_of_attack(theta, x_velocity, y_velocity):
    return np.arctan2(x_velocity, y_velocity) - theta

def dynamicViscosity(altitude):
    return 1.458*10**(-6)*((temperature(altitude)+273.15)**(3/2))/(temperature(altitude)+273.15+110.4)

def kinematicViscosity(altitude):
    return dynamicViscosity(altitude)/airDensity(altitude)

def reynoldsNumber(windSpeed, rocketVelocityX, rocketVelocityY, altitude, characteristicLength):
    return airSpeed(windSpeed, rocketVelocityX, rocketVelocityY) * characteristicLength / kinematicViscosity(altitude)

def zeroAngleDragCoefficient(time, windSpeed, rocketVelocityX, rocketVelocityY, altitude, L):
    return C_D_friction(windSpeed, rocketVelocityX, rocketVelocityY, altitude, L) + nosePressureDrag(windSpeed, rocketVelocityX, rocketVelocityY, altitude) + finPressureDrag(time, windSpeed, rocketVelocityX, rocketVelocityY, altitude) + baseDragCoefficient(time, windSpeed, rocketVelocityX, rocketVelocityY, altitude) + parasiticDrag(windSpeed, rocketVelocityX, rocketVelocityY, altitude)

def C_D_friction(windSpeed, rocketVelocityX, rocketVelocityY, altitude, L):
    R = reynoldsNumber(windSpeed, rocketVelocityX, rocketVelocityY, altitude, L)
    M = machNumber(windSpeed, rocketVelocityX, rocketVelocityY, altitude)
    R_crit = 51*(surfaceRoughness/L)**-1.039
    if R < 10**4:
        C_f = 1.48*10**-2
    elif R < R_crit:
        C_f = 1/(1.5*np.log(R)-5.6)**2
    else:
        C_f = 0.032*(surfaceRoughness/L)**0.2

    C_Mach = (1-0.1*M**2)
    C_f_rough = 0.032*(surfaceRoughness/L)**0.2 * C_Mach
    C_f_component = max(C_f, C_f_rough)
    K_body = 1 + 1/(2*finenessRatio)
    K_fin = 1 + 2*finThickness/meanAerodynamicChordLengthOfFin
    C_D_body = C_f_component * K_body * (A_wet_body/A_ref)
    C_D_fins = C_f_component * K_fin * (A_wet_fins/A_ref)
    return C_D_body + C_D_fins

def nosePressureDrag(windSpeed, rocketVelocityX, rocketVelocityY, altitude):
    kappa=(1/1) #rho_t/rho in current case both are equal so for performance, it will be simplified.
    M = machNumber(windSpeed, rocketVelocityX, rocketVelocityY, altitude)
    epsilon = np.arctan(bodyDiameter/(2*noseLength))
    gamma = 1.4
    C_D_At_M1 = np.sin(epsilon)
    slope = 4/(gamma+1) * (1-0.5*C_D_At_M1)
    if M < 0.8:
        coneDragCoefficient = 0.8*np.sin(jointAngle)**2
    elif M < 1.2:
        coneDragCoefficient = (3*slope+C_D_At_M1-2*np.sin(jointAngle)**2)*(M-0.8)+0.8*np.sin(jointAngle)**2
    else:
        coneDragCoefficient = slope*M + C_D_At_M1
    return  (0.72*(kappa-0.5)**2 + 0.82)* coneDragCoefficient #the correction factor is only used for ogival shapes, which we are, therefor for now well assume it  

def finPressureDrag(time, windSpeed, rocketVelocityX, rocketVelocityY, altitude):
    M = machNumber(windSpeed, rocketVelocityX, rocketVelocityY, altitude)
    if M < 0.9: #rounded edge
        perpindicularLeadingEdgeDrag = (1-M**2)**(-0.417) -1
    elif M < 1:
        perpindicularLeadingEdgeDrag = 1 - 1.785*(M-0.9)
    else:
        perpindicularLeadingEdgeDrag = 1.214 - (0.502/M**2) + (0.1095/M**4)
    leadingEdgeDrag = perpindicularLeadingEdgeDrag * np.cos(leadingEdgeAngle)**2 #angled fin, sweep angle
    trailingEdgeDrag = (1/2) * baseDragCoefficient(time, windSpeed, rocketVelocityX, rocketVelocityY, altitude) #rounded edge
    C_D_fin = leadingEdgeDrag + trailingEdgeDrag
    A_FIN = numberOfFins * finThickness * finSpan
    return (A_FIN/A_ref) * C_D_fin

def baseDragCoefficient(time, windSpeed, rocketVelocityX, rocketVelocityY, altitude):
    if time < timeCurve[-1]:
        A_motor = (np.pi/4) * motorDiameter ** 2
    else:
        A_motor = 0
    M = machNumber(windSpeed, rocketVelocityX, rocketVelocityY, altitude)
    if M < 1:
        CD_base = 0.12 + 0.13*M**2
    else:
        CD_base = 0.25/M
    return CD_base * (A_ref-A_motor)/A_ref

def parasiticDrag(windSpeed, rocketVelocityX, rocketVelocityY, altitude):
    M = machNumber(windSpeed, rocketVelocityX, rocketVelocityY, altitude)
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
    f = 0.25/((np.log10((launchGuideRoughness/(3.7*launchGuideInnerDiameter)) + (5.74/(reynoldsNumber(windSpeed, rocketVelocityX, rocketVelocityY, altitude, launchGuideInnerDiameter)**0.9))))**2)
    C_D_tube = f * (launchGuideLength/launchGuideInnerDiameter)
    return (C_D_tube * inner_Area + 0.7*(C_D_parasitic + C_D_base)*A_parasitic)/A_ref

def axialDragCoefficient(AOA):
    if abs(AOA) <= 17:
        scalefunction = (-3/24565)*abs(AOA)**3 + (9/2890)*AOA**2 + 1
    else:
        scalefunction = (13/1945085)*(AOA-90)**3 + (2847/3890170)*(AOA-90)**2
    return zeroAngleDragCoefficient() * scalefunction

def axialDragForce():
    rho = airDensity()
    V = airSpeed()
    C_A = axialDragCoefficient
    return (1/2)*rho*V**2 * A_ref * C_A

def normalForceFinCoefficientDerivative(AOA, windSpeed, rocketVelocityX, rocketVelocityY, altitude):
    M = machNumber(windSpeed, rocketVelocityX, rocketVelocityY, altitude)
    midchordSweepAngle = np.arctan(finSpan/(sweepLength+tipChord/2-rootChord/2))
    beta = np.sqrt(abs(M**2-1))
    C_Nalpha0 = 2*np.pi/beta
    AR = 2*(finSpan**2)/A_fin
    F_D = AR/((1/(2*np.pi))*C_Nalpha0*np.cos(midchordSweepAngle))
    C_N_1fin = (C_Nalpha0 * F_D * (A_fin/A_ref) * np.cos(midchordSweepAngle))/(2+F_D*np.sqrt(1+4/(F_D**2)))    #(2*np.pi*((finSpan**2)/A_ref))/(2+((beta*finSpan**2)/(A_fin*np.cos(midchordSweepAngle))))
    C_N_fins = (numberOfFins/2) * C_N_1fin
    K_TB = 1 + (bodyRadius/(finSpan+bodyRadius))
    return K_TB * C_N_fins

def windSpeed(altitude):
    return initWindSpeed * (altitude/initialAltitude)


def normalForceCoefficientDerivative(AOA, windSpeed, rocketVelocityX, rocketVelocityY, altitude):
    C_N_finsWithInterference = normalForceFinCoefficientDerivative(AOA, windSpeed, rocketVelocityX, rocketVelocityY, altitude)
    C_N_nose = 2*np.cos(AOA)
    return C_N_nose + C_N_finsWithInterference

def pitchMomentCoefficientDerivative(AOA, windSpeed, rocketVelocityX, rocketVelocityY, altitude):
    C_m_nose = (2*np.sin(AOA))/(A_ref*bodyDiameter)*(noseLength*A_ref-noseVolume)
    M = machNumber(windSpeed, rocketVelocityX, rocketVelocityY, altitude)
    AR = 2*finSpan**2/A_fin
    beta = np.sqrt(abs(M**2-1))

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
    C_m_fins = normalForceFinCoefficientDerivative(AOA, windSpeed, rocketVelocityX, rocketVelocityY, altitude) * x_fin/bodyDiameter
    return C_m_nose+C_m_fins

def normalForce(AOA, windSpeed, rocketVelocityX, rocketVelocityY, altitude):
    rho = airDensity(altitude)
    V = airSpeed(windSpeed, rocketVelocityX, rocketVelocityY)
    C_N_alpha = normalForceCoefficientDerivative(AOA,windSpeed, rocketVelocityX, rocketVelocityY, altitude)
    return (1/2)*rho*V**2 * A_ref * bodyDiameter * C_N_alpha * AOA

def pitchMoment(AOA,windSpeed, rocketVelocityX, rocketVelocityY, altitude):
    rho = airDensity(altitude)
    V = airSpeed(windSpeed,rocketVelocityX,rocketVelocityY)
    C_m_alpha = pitchMomentCoefficientDerivative(AOA, windSpeed, rocketVelocityX, rocketVelocityY, altitude)
    return (1/2)*rho*V**2 * A_ref * C_m_alpha * AOA

def centerOfPressure(AOA,windSpeed, rocketVelocityX, rocketVelocityY, altitude):
    if AOA != 0:
        return (pitchMomentCoefficientDerivative(AOA,windSpeed, rocketVelocityX, rocketVelocityY, altitude)/normalForceCoefficientDerivative(AOA,windSpeed, rocketVelocityX, rocketVelocityY, altitude))*bodyDiameter
    else:
        return CP

def windSpeedatAltitude(altitude):
    return initWindSpeed + 0*altitude

def moment(time, theta, rocketVelocityX, rocketVelocityY, altitude):
    AOA = angle_of_attack(theta, rocketVelocityX, rocketVelocityY)
    windSpeed = windSpeedatAltitude(altitude)
    CG = centerOfGravity(time)
    return normalForce(AOA, windSpeed, rocketVelocityX, rocketVelocityY, altitude)*(centerOfPressure(AOA, windSpeed, rocketVelocityX, rocketVelocityY, altitude)-CG)

# #Functions






def moment_of_inertia(current_mass, time):
    I_center = (1/12) * current_mass * rocketLength ** 2 #assuming uniform density for now
    d = centerOfGravity(time) - rocketLength/2 #distance between center and CG
    return I_center + current_mass * d**2 #parallel axis theorem

def normal_force(theta, y_pos, x_velocity, y_velocity):
    dynamicPressure = 0.5*air_density(y_pos)*netSpeed(x_velocity, y_velocity)**2 #Assuming incompressible air for now
    Area = np.pi * radius**2
    return dynamicPressure * Area * CN_alpha * angle_of_attack(theta, x_velocity,y_velocity) 

def yaw_torque(theta, y_pos, x_velocity, y_velocity, time):
    return normal_force(theta, y_pos, x_velocity, y_velocity) * (CP-centerOfGravity(time))

def alpha(time, theta, y_pos, x_velocity, y_velocity):
    return moment(time, theta, x_velocity, y_velocity, y_pos) / moment_of_inertia(mass(time),time)

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

def x_component(magnitude, theta):
    return magnitude * np.sin(theta)

def y_component(magnitude, theta):
    return magnitude * np.cos(theta)

def acceleration_x(time, y_pos, x_velocity, y_velocity, theta):
    return netForce_x(time, y_pos, x_velocity, y_velocity, theta) / mass(time)

def acceleration_y(time, y_pos, x_velocity, y_velocity, theta):
    return netForce_y(time, y_pos, x_velocity, y_velocity, theta) / mass(time)

def netForce_x(time, y_pos, x_velocity, y_velocity, theta):
    if y_pos == 0:
        return x_component(thrust_magnitude(time), theta)
    else:
        return x_component(thrust_magnitude(time), theta) - x_component(drag_magnitude(y_pos, x_velocity, y_velocity), np.arctan2(x_velocity, y_velocity))    

def netForce_y(time, y_pos, x_velocity, y_velocity, theta):
    if y_pos == 0:
        return y_component(thrust_magnitude(time), theta)
    else:
        return y_component(thrust_magnitude(time), theta) - weight(time) - y_component(drag_magnitude(y_pos, x_velocity, y_velocity),np.arctan2(x_velocity, y_velocity))

def drag_magnitude(y_pos, x_velocity, y_velocity):
    speed = netSpeed(x_velocity, y_velocity)
    if y_velocity >= 0:
        A, cd = 3.14 * radius ** 2, dragCoefficient        
    else:
        A, cd = parachuteArea, parachuteDragCoefficient
    return (0.5 * air_density(y_pos) * cd * A * speed**2)

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

def centerOfGravity(time):
    return (CG_dry*dryMass+CT*(mass(time)-dryMass))/(mass(time))
    

def rk4_step(time, x_pos, y_pos, x_velocity, y_velocity, theta, omega):
    k1_x = x_velocity
    k1_y = y_velocity
    k1_vx = acceleration_x(time, y_pos, x_velocity, y_velocity, theta)
    k1_vy = acceleration_y(time, y_pos, x_velocity, y_velocity, theta)
    k1_0 = omega
    k1_w = alpha(time, theta, y_pos, x_velocity, y_velocity)

    k2_x = x_velocity + k1_vx * timeStep/2
    k2_y = y_velocity + k1_vy * timeStep/2
    k2_vx = acceleration_x(time + timeStep/2, y_pos + k1_y * timeStep/2, x_velocity + k1_vx * timeStep/2, y_velocity + k1_vy * timeStep/2, theta + k1_0 * timeStep/2)
    k2_vy = acceleration_y(time + timeStep/2, y_pos + k1_y * timeStep/2, x_velocity + k1_vx * timeStep/2, y_velocity + k1_vy * timeStep/2, theta + k1_0 * timeStep/2)
    k2_0 = omega + k1_w * timeStep/2
    k2_w = alpha(time + timeStep/2, theta + k1_0 * timeStep/2, y_pos + k1_y * timeStep/2, x_velocity + k1_vx * timeStep/2, y_velocity + k1_vy * timeStep/2)

    k3_x = x_velocity + k2_vx * timeStep/2
    k3_y = y_velocity + k2_vy * timeStep/2
    k3_vx = acceleration_x(time + timeStep/2, y_pos + k2_y * timeStep/2, x_velocity + k2_vx * timeStep/2, y_velocity + k2_vy * timeStep/2, theta + k2_0 * timeStep/2)
    k3_vy = acceleration_y(time + timeStep/2, y_pos + k2_y * timeStep/2, x_velocity + k2_vx * timeStep/2, y_velocity + k2_vy * timeStep/2, theta + k2_0 * timeStep/2)
    k3_0 = omega + k2_w * timeStep/2
    k3_w = alpha(time + timeStep/2, theta + k2_0 * timeStep/2, y_pos + k2_y * timeStep/2, x_velocity + k2_vx * timeStep/2, y_velocity + k2_vy * timeStep/2)

    k4_x = x_velocity + k3_vx * timeStep
    k4_y = y_velocity + k3_vy * timeStep
    k4_vx = acceleration_x(time + timeStep, y_pos + k3_y * timeStep, x_velocity + k3_vx * timeStep, y_velocity + k3_vy * timeStep, theta + k3_0 * timeStep)
    k4_vy = acceleration_y(time + timeStep, y_pos + k3_y * timeStep, x_velocity + k3_vx * timeStep, y_velocity + k3_vy * timeStep, theta + k3_0 * timeStep)
    k4_0 = omega + k3_w * timeStep
    k4_w = alpha(time + timeStep, theta + k3_0 * timeStep, y_pos + k3_y * timeStep, x_velocity + k3_vx * timeStep, y_velocity + k3_vy * timeStep)


    new_x_pos = x_pos + (timeStep/6) * (k1_x + 2*k2_x + 2*k3_x + k4_x)
    new_y_pos = y_pos + (timeStep/6) * (k1_y + 2*k2_y + 2*k3_y + k4_y)
    new_x_velocity = x_velocity + (timeStep/6) * (k1_vx + 2*k2_vx + 2*k3_vx + k4_vx)
    new_y_velocity = y_velocity + (timeStep/6) * (k1_vy + 2*k2_vy + 2*k3_vy + k4_vy)
    new_theta = theta + (timeStep/6) * (k1_0 + 2*k2_0 + 2*k3_0 + k4_0)
    new_omega = omega + (timeStep/6) * (k1_w + 2*k2_w + 2*k3_w + k4_w)

    if(on_rod(new_x_pos, new_y_pos)):
        new_theta = theta
        new_omega = 0
    
    if new_y_velocity < 0:
        new_theta = theta
        new_omega = 0

    new_theta = (new_theta + np.pi) % (2*np.pi) - np.pi

    return new_x_pos, new_y_pos, new_x_velocity, new_y_velocity, new_theta, new_omega
    
def run_simulation():
    global time, x_pos, y_pos, x_velocity, y_velocity, theta, omega
    while y_pos >= -0.01:

        current_mass = mass(time)
        current_thust = thrust_magnitude(time)
        current_thrust_x = x_component(current_thust, theta)
        current_thrust_y = y_component(current_thust, theta)
        current_weight = current_mass * gravity
        current_drag = drag_magnitude(y_pos, x_velocity, y_velocity)
        current_drag_x = x_component(current_drag, np.arctan2(x_velocity, y_velocity))
        current_drag_y = y_component(current_drag, np.arctan2(x_velocity, y_velocity))
        current_normalForce = normal_force(theta, y_pos, x_velocity, y_velocity)
        current_inertia = moment_of_inertia(current_mass, centerOfGravity(time))
        
        if y_pos == 0:
            current_netForce_x = current_thrust_x
            current_netForce_y = current_thrust_y 
        else:
            current_netForce_x = current_thrust_x - current_drag_x
            current_netForce_y = current_thrust_y - current_weight - current_drag_y

        if on_rod(x_pos, y_pos):
            current_torque = 0
        else:
            current_torque = yaw_torque(theta,y_pos,x_velocity,y_velocity,time)

        current_acceleration_x = current_netForce_x / current_mass
        current_acceleration_y = current_netForce_y / current_mass
        current_alpha = current_torque/current_inertia

        timeList.append(time)
        x_posList.append(x_pos)
        y_posList.append(y_pos)
        x_velocityList.append(x_velocity)  
        y_velocityList.append(y_velocity)
        thetaList.append(np.degrees(theta))
        omegaList.append(omega)
        alphaList.append(current_alpha)
        torqueList.append(current_torque)
        inertiaList.append(current_inertia)
        normalForceList.append(current_normalForce)
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
        new_x_pos, new_y_pos, new_x_velocity, new_y_velocity, new_theta, new_omega = rk4_step(time, x_pos, y_pos, x_velocity, y_velocity, theta, omega)
        x_pos = new_x_pos
        y_pos = new_y_pos
        x_velocity = new_x_velocity
        y_velocity = new_y_velocity
        theta = new_theta
        omega = new_omega
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
        print(f"Time: {timeList[i]:.3f} s, \
              Mass: {massList[i]:.3f} kg, \
              x_pos: {x_posList[i]:.3f} m, \
              Altitude: {y_posList[i]:.3f} m, \
              x_vel: {x_velocityList[i]:.3f} m/s, \
              Velocity: {y_velocityList[i]:.3f} m/s, \
              x_acc: {x_accelerationList[i]:.3f} m/s^2 \
              acceleration: {y_accelerationList[i]:.3f} m/s^2, \
              x_net: {x_netForceList[i]:.3f} N \
              Net Force: {y_netForceList[i]:.3f} N, \
              Weight: {weightList[i]:.3f} N, \
              x_thrust: {x_thrustList[i]:.3f} N \
              Thrust: {y_thrustList[i]:.3f} N, \
              x_drag: {x_dragList[i]:.3f} N, \
              Drag: {y_dragList[i]:.3f} N, \
              theta: {thetaList[i]:.3f} degrees, \
              omega: {omegaList[i]:.3f} rad/s, \
              alpha: {alphaList[i]:.3f} rad/s^2, \
              torque: {torqueList[i]:.3f} N*m, \
              inertia: {inertiaList[i]:.3f} kg*m^2, \
              normalF: {normalForceList[i]:.3f} N, \
              Air Density: {densityList[i]:.3f} kg/m^3")
    rodVelocity = 0
    for i in range(len(y_posList)-1):
        dist1=np.sqrt(y_posList[i]**2+x_posList[i]**2)
        dist2=np.sqrt(y_posList[i+1]**2+x_posList[i+1]**2)
        if dist1 < rodLength <= dist2:
            rodVelocity = np.sqrt((y_velocityList[i]**2 + 2 * y_accelerationList[i] * (rodLength*np.cos(launchAngle) - y_posList[i]))+(x_velocityList[i]**2 + 2 * x_accelerationList[i] * (rodLength*np.sin(launchAngle) - x_posList[i])))
            break
    print(f"Velocity off rod: {rodVelocity:.3f} m/s, Apogee: {max(y_posList):.3f} m, Max Velocity: {max(y_velocityList):.3f} m/s, Max acceleration: {max(y_accelerationList):.3f} m/s^2, Time to Apogee: {timeList[y_posList.index(max(y_posList))]:.3f} s, Flight Time: {timeList[-1]:.3f} s, Ground hit velocity: {y_velocityList[-1]:.3f} m/s, Range from launch: {x_posList[-1]:.3f} m")
    sys.stdout.flush()
    plot_results()