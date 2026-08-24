import numpy as np
import random

#verityables
initialTemperature = 20 #degrees Celsius
initialAltitude = 220 #meters above sea level
initialPressure = 101550 #Pa
airMolarMass = 0.0289644 #kg/mol
universeGasConstant = 8.31447 #J/(mol*K)
SpecificGasConstant = 287.05 #J/(kg*K)
temperatureLapseRate = 0.0065 #K/m
gravity = 9.806 #m/s^2
AdiabaticIndex = 1.4 #unitless & constant for our purpose, value for air
windSpeed = 0 #m/s
windheading = 0 #radians, 0 is wind to the north, pi/2 is wind to the east, pi is wind to the south, 3pi/2 is wind to the west
surfaceRoughness = 60*10**-9 #meters, first number is nanometer, you have to look in a table or guesstimate for this value. 
numberOfFins = 3
rocketLength = 0.4 #meters
finThickness = 0.002 #meters
rootChord = 0.0508 #meters
tipChord = 0.0508 #meters
sweepLength = 0.0254 #meters
finSpan = 0.03 #meters
leadingEdgeToNosecone = 0.349 #meters 
bodyDiameter = 0.025 #meters
bodyRadius = bodyDiameter/2 #meters
motorDiameter = 0.018 # meters also just the inner body tube innerdiameter
noseLength = 0.1 #meters
A_fin = numberOfFins * finThickness * finSpan
finenessRatio = rocketLength/bodyDiameter #unitless, ratio of the length of the rocket to the diameter of the rocket
meanAerodynamicChordLengthOfFin = (2/3) * (rootChord + tipChord - (rootChord*tipChord)/(rootChord + tipChord)) #meters, this is a typical value for a fin
A_wet_nose = (np.pi/(4*bodyRadius**2))*(2*bodyRadius*noseLength*(bodyRadius**2 - noseLength**2)+((bodyRadius**2 + noseLength**2)**2)*np.arcsin(bodyRadius/noseLength)) #meters^2, this is the wetted area of an ogive nosecone of the rocket
noseVolume = (np.pi/(24*bodyRadius**3))*(6*noseLength*bodyRadius**5+6*bodyRadius*noseLength**5+4*bodyRadius**3*noseLength**3+6*(bodyRadius**2+noseLength**2)*(bodyRadius**2+noseLength**2)**2*np.arctan(bodyRadius/noseLength))
A_wet_body = np.pi * bodyDiameter * rocketLength + A_wet_nose #meters^2, this is the wetted area of the body of the rocket
A_wet_fins = numberOfFins*((rootChord+tipChord)+finThickness*(rootChord+tipChord+np.sqrt(sweepLength**2 + finSpan**2)+np.sqrt((sweepLength+tipChord-rootChord)**2+finSpan**2))) #meters^2, this is the wetted area of the fins of the rocket
A_ref = np.pi * (bodyRadius**2) #meters^2, this is the reference area of the rocket, which is the cross-sectional area of the rocket body
jointAngle = 0 #radians, this is the angle of the joint between the body and the fin, which is typically 0 for a rocket with fins that are perpendicular to the body
launchGuideLength = 0.035 #meters
launchGuideOuterDiameter = 0.007 # meters
launchGuideInnerDiameter = 0.0005 #meters

def temperature(altitude):
    if altitude+initialAltitude < 11000:
        return initialTemperature - temperatureLapseRate * (altitude)
    else:
        return initialTemperature - temperatureLapseRate * (11000-initialAltitude)

def pressure(altitude):
    return initialPressure * (temperature(altitude)+273.15)/(initialTemperature+273.15)**(gravity*airMolarMass/(universeGasConstant*temperatureLapseRate))

def dynamicPressure(windSpeed, windHeading, rocketVelocityX, rocketVelocityY, altitude):
    return 0.5 * airDensity(altitude) * airSpeed(windSpeed, windHeading, rocketVelocityX, rocketVelocityY)**2

def airDensity(altitude):
    return pressure(altitude)/(SpecificGasConstant*(temperature(altitude)+273.15))

def speedOfSound(altitude):
    return np.sqrt(AdiabaticIndex*SpecificGasConstant*(temperature(altitude)+273.15))

def machNumber(windSpeed, windHeading, rocketVelocityX, rocketVelocityY, altitude):
    return airSpeed(windSpeed, windHeading, rocketVelocityX, rocketVelocityY)/speedOfSound(altitude)

def airSpeed(windspeed, windHeading, rocketVelocityX, rocketVelocityY):
    return np.sqrt((rocketVelocityX-windspeed*np.sin(windHeading))**2 + (rocketVelocityY**2) + (windspeed*np.cos(windHeading))**2)

def rocketSpeed(rocketVelocityX, rocketVelocityY):
    return np.sqrt(rocketVelocityX**2 + rocketVelocityY**2)

def angle_of_attack(theta, x_velocity, y_velocity):
    return np.arctan2(x_velocity, y_velocity) - theta

def dynamicViscosity(altitude):
    return 1.458*10**(-6)*((temperature(altitude)+273.15)**(3/2))/(temperature(altitude)+273.15+110.4)

def kinematicViscosity(altitude):
    return dynamicViscosity(altitude)/airDensity(altitude)

def reynoldsNumber(windSpeed, windHeading, rocketVelocityX, rocketVelocityY, altitude, characteristicLength):
    return airSpeed(windSpeed, windHeading, rocketVelocityX, rocketVelocityY) * characteristicLength / kinematicViscosity(altitude)

def zeroAngleDragCoefficient():
    return C_D_friction() + nosePressureDrag() + finPressureDrag() + baseDragCoefficient() + parasiticDrag()

def C_D_friction(windSpeed, windHeading, rocketVelocityX, rocketVelocityY, altitude, L):
    R = reynoldsNumber(windSpeed, windHeading, rocketVelocityX, rocketVelocityY, altitude, L)
    R_crit = 51*(surfaceRoughness/L)**-1.039
    if R < 10**4:
        C_f = 1.48*10**-2
    elif R < R_crit:
        C_f = 1/(1.5*np.log(R)-5.6)**2
    else:
        C_f = 0.032*(surfaceRoughness/L)**0.2
    C_f_c = C_f(1-0.1*(machNumber(rocketVelocityX, rocketVelocityY, altitude)**2))
    return C_f_c*(1+(1/(2*finenessRatio)) * A_wet_body + (1+(2*finThickness)/meanAerodynamicChordLengthOfFin) * A_wet_fins)/(A_ref)

def nosePressureDrag(windSpeed, windHeading, rocketVelocityX, rocketVelocityY, altitude):
    kappa=(1/1) #rho_t/rho in current case both are equal so for performance, it will be simplified.
    M = machNumber(windSpeed, windHeading, rocketVelocityX, rocketVelocityY, altitude)
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

def finPressureDrag(windSpeed, windHeading, rocketVelocityX, rocketVelocityY, altitude):
    M = machNumber(windSpeed, windHeading, rocketVelocityX, rocketVelocityY, altitude)
    if M < 0.9: #rounded edge
        perpindicularLeadingEdgeDrag = (1-M**2)**(-0.417) -1
    elif M < 1:
        perpindicularLeadingEdgeDrag = 1 - 1.785*(M-0.9)
    else:
        perpindicularLeadingEdgeDrag = 1.214 - (0.502/M**2) + (0.1095/M**4)
    leadingEdgeAngle = np.arctan(finSpan/sweepLength)
    leadingEdgeDrag = perpindicularLeadingEdgeDrag * np.cos(leadingEdgeAngle)**2 #angled fin, sweep angle
    trailingEdgeDrag = (1/2)* baseDragCoefficient() #rounded edge
    C_D_fin = leadingEdgeDrag + trailingEdgeDrag
    return (A_fin/A_ref) * C_D_fin

def baseDragCoefficient(time, timecurve):
    if time < timecurve[-1]:
        A_motor = (np.pi/4) * motorDiameter ** 2
    else:
        A_motor = 0
    M = machNumber()
    if M < 1:
        CD_base = 0.12 + 0.13*M**2
    else:
        CD_base = 0.25/M
    return CD_base * (A_ref-A_motor)/A_ref

def parasiticDrag():
    M = machNumber()
    if M < 1:
        qRatio = 1 + (M**2/4) + (M**4/40)
    else:
        qRatio = 1.84 - (0.76/M**2) + (0.166/M**4) + (0.035/M**6)
    stagDrag = 0.85 * qRatio
    C_D_parasitic = max(1.3-0.3*(launchGuideLength/launchGuideOuterDiameter), 1) * stagDrag
    A_parasitic = (np.pi/4) * (launchGuideOuterDiameter**2-launchGuideInnerDiameter**2 * max(1-launchGuideLength/launchGuideOuterDiameter,0))
    return (A_parasitic/A_ref) * C_D_parasitic

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

def normalForceFinCoefficient(AOA):
    M = machNumber()
    midchordSweepAngle = np.pi - np.arctan(finSpan/(sweepLength+tipChord/2-rootChord/2))
    beta = np.sqrt(abs(M**2-1))
    C_N_1fin = (2*np.pi*((finSpan**2)/A_ref))/(1+np.sqrt(1+((beta*finSpan**2)/(A_fin*np.cos(midchordSweepAngle)))**2))
    C_N_fins = (numberOfFins/2) * C_N_1fin
    K_TB = 1 + bodyRadius/(finSpan+bodyRadius)
    return K_TB * C_N_fins * np.sin(AOA)


def normalForceCoefficient(AOA):
    C_N_finsWithInterference = normalForceFinCoefficient()
    C_N_nose = 2*np.sin(AOA)
    return C_N_nose + C_N_finsWithInterference

def pitchMomentCoefficient(AOA):
    C_m_nose = (2*np.sin(AOA))/(A_ref*bodyDiameter)*(noseLength*A_ref-noseVolume)
    C_m_body = (2*np.sin(AOA))/(A_ref*bodyDiameter)*(rocketLength-noseLength)*(A_ref - np.pi*bodyRadius**2)
    M = machNumber()
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

    coefficients = np.linealg.solve(A,b)

    if M <= 0.5:
        x_f = (sweepLength/3)*((rootChord+2*tipChord)/(rootChord+tipChord))+(1/6)*((rootChord**2+tipChord**2+rootChord*tipChord)/(rootChord+tipChord))
    elif M <=2:
        x_f = coefficients[0] + coefficients[1]*M + coefficients[2]*M**2 + coefficients[3]*M**3 + coefficients[4]*M**4 + coefficients[5]*M**5
    else:
        x_f = meanAerodynamicChordLengthOfFin*((AR*beta - 0.67)/(2*AR*beta-1))
    x_fin = leadingEdgeToNosecone + x_f
    C_m_fins = normalForceFinCoefficient() * x_fin/bodyDiameter
    return C_m_nose+C_m_body+C_m_fins

def normalForce(AOA):
    rho = airDensity()
    V = airSpeed()
    C_m = pitchMomentCoefficient(AOA)
    return (1/2)*rho*V**2 * A_ref * bodyDiameter * C_m

def pitchMoment(AOA):
    rho = airDensity()
    V = airSpeed()
    C_N = normalForceCoefficient(AOA)
    return (1/2)*rho*V**2 * A_ref * C_N

def centerOfPressure(AOA):
    return (pitchMomentCoefficient(AOA)/normalForceCoefficient(AOA))*bodyDiameter

def moment(centerOfGravity,AOA):
    return normalForce(AOA)*(centerOfPressure(AOA)-centerOfGravity)

print(normalForceCoefficient(0))
print(moment(0.26,0))