#sim.py - V2.3

import matplotlib.pyplot as plt 
import numpy as np
import csv
import sys
import io

class RocketSimNew:

    def import_Motor_Data(self, filename):
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

    #Constants
    def __init__(self, motor_source="Motors/Estes_C6.csv", **overrides):
        self.gravity = 9.81 # m/s^2
        self.timeStep = 0.05 # seconds
        self.initialTemperature = 20 #degrees Celsius
        self.initialAltitude = 220 #meters above sea level
        self.initialPressure = 101550 #Pa
        self.airMolarMass = 0.0289644 #kg/mol
        self.universeGasConstant = 8.31447 #J/(mol*K)
        self.SpecificGasConstant = 287.05 #J/(kg*K)
        self.temperatureLapseRate = 0.0065 #K/m
        self.hellMannCoefficient = 0.34 #constant for Neutral air above human-inhibited area
        self.AdiabaticIndex = 1.4 #unitless & constant for our purpose, value for air
        self.initWindSpeed = 0 #m/s
        self.windHeading = 0 #radians, 0 is wind to the north, pi/2 is wind to the east, pi is wind to the south, 3pi/2 is wind to the west
        self.surfaceRoughness = 60*10**-6 #meters, first number is micrometer, you have to look in a table or guesstimate for this value. 
        self.launchGuideRoughness = 60*10**-9
        self.numberOfFins = 3
        self.rocketLength = 0.425 #meters
        self.finThickness = 0.002 #meters
        self.rootChord = 0.0508 #meters
        self.tipChord = 0.0508 #meters
        self.sweepLength = 0.0254 #meters
        self.leadingEdgeAngle = np.radians(40.3)
        self.finSpan = 0.03 #meters
        self.Y1 = np.sqrt(self.finSpan**2+(self.sweepLength+self.tipChord/2-self.rootChord/2)**2) #supposed to be lambda or the distance on Fin at the mid-chord lines for set 1 (m)
        self.leadingEdgeToNosecone = 0.349 #meters 
        self.bodyDiameter = 0.025 #meters
        self.bodyRadius = self.bodyDiameter/2 #meters
        self.motorDiameter = 0.018 # meters also just the inner body tube innerdiameter
        self.noseLength = 0.1 #meters
        self.A_fin = (self.finSpan/2)*(self.tipChord+self.rootChord)
        self.finenessRatio = self.rocketLength/self.bodyDiameter #unitless, ratio of the length of the rocket to the diameter of the rocket
        self.meanAerodynamicChordLengthOfFin = (2/3) * (self.rootChord + self.tipChord - (self.rootChord*self.tipChord)/(self.rootChord + self.tipChord)) #meters, this is a typical value for a fin
        self.A_wet_nose = (np.pi/(2*self.bodyRadius**2))*(self.bodyRadius*self.noseLength*(self.bodyRadius**2-self.noseLength**2)+((self.bodyRadius**2+self.noseLength**2)**2)*np.atan(self.bodyRadius/self.noseLength)) #meters^2, this is the wetted area of an ogive nosecone of the rocket
        self.noseVolume = (np.pi/(24*self.bodyRadius**3))*(6*self.noseLength*self.bodyRadius**5+6*self.bodyRadius*self.noseLength**5+4*(self.bodyRadius**3)*(self.noseLength**3)+6*(self.bodyRadius**2-self.noseLength**2)*(self.bodyRadius**2+self.noseLength**2)**2*np.arctan(self.bodyRadius/self.noseLength))
        self.A_wet_body = np.pi * self.bodyDiameter * (self.rocketLength-self.noseLength) + self.A_wet_nose #meters^2, this is the wetted area of the body of the rocket
        self.A_wet_fins = self.numberOfFins*(self.finSpan*(self.rootChord+self.tipChord)+self.finThickness*(self.rootChord+self.tipChord+np.sqrt(self.sweepLength**2 + self.finSpan**2)+np.sqrt((self.sweepLength+self.tipChord-self.rootChord)**2+self.finSpan**2))) #meters^2, this is the wetted area of the fins of the rocket
        self.A_ref = np.pi * (self.bodyRadius**2) #meters^2, this is the reference area of the rocket, which is the cross-sectional area of the rocket body
        self.jointAngle = 0 #radians, this is the angle of the joint between the body and the fin, which is typically 0 for a rocket with fins that are perpendicular to the body
        self.launchGuideLength = 0.035 #meters
        self.launchGuideOuterDiameter = 0.007 # meters
        self.launchGuideInnerDiameter = 0.005 #meters

        self.finSets = 1
        self.transitions = 0
        self.noseType = "ogive" #cone, ogive, paraboloid, ellipsoid


        #EDITABLE
        self.dryMass = 0.0605 # kg
        self.fuelMass = 0.011 # kg
        self.radius = 0.0127 # m
        self.dragCoefficient = 0.634 # dimensionless
        self.parachuteArea = 0.0707 # m^2
        self.parachuteDragCoefficient = 0.80 # dimensionless
        self.rodLength = 1.0 # m
        self.launchAngle = np.radians(5) #radians of x degrees
        self.windSpeed = 0 # m/s
        self.rocketLength = 0.425 # m
        self.CG_dry = 0.24 # m from the nose tip
        self.CG_wet = 0.26
        self.CP = self.initcenterOfPressure() # m from the nose tip
        self.CT = (self.CG_wet*(self.dryMass+self.fuelMass)-self.CG_dry*self.dryMass)/(self.fuelMass)

        
    
        self.CN_alpha = 11.97 #pulled from OR
    
        #fileName = input("Enter the CSV file containing motor data: ") # CSV file containing motor data
    
        #motor_csv = import_Motor_Data(fileName)
        motor_csv = self.import_Motor_Data("Motors/Estes_C6.csv")
        self.timeCurve = motor_csv[0]
        self.thrustCurve = motor_csv[1]
    
        #Init Variables
        self.time = 0
        self.y_pos = 0
        self.y_velocity = 0
        self.x_pos = 0
        self.x_velocity = 0
        self.theta = self.launchAngle
        self.omega = 0
        self.leftRod = False
        for key, value in overrides.items():
            setattr(self, key, value)
        #Init Lists
        self.timeList = []
        self.massList = []
        self.x_posList = []
        self.y_posList = []
        self.x_velocityList = []
        self.y_velocityList = []
        self.thetaList = []
        self.omegaList = []
        self.alphaList = []
        self.torqueList = []
        self.inertiaList = []
        self.normalForceList = []
        self.x_accelerationList = []
        self.y_accelerationList = []
        self.x_netForceList = []
        self.y_netForceList = []
        self.weightList = []
        self.x_thrustList = []
        self.y_thrustList = []
        self.x_dragList = []
        self.y_dragList = []
        self.densityList = []

    def initcenterOfPressure(self):
        CN_n = 2
        CN_f1 = (1+(self.bodyRadius)/(self.finSpan+self.bodyRadius)) * ((4*self.numberOfFins*(self.finSpan/self.bodyDiameter)**2)/(1+np.sqrt(1+((2*self.Y1)/(self.rootChord+self.tipChord))**2)))
        Pf1 = CN_f1 * (self.leadingEdgeToNosecone+(self.sweepLength*(self.rootChord+2*self.tipChord)/(3*(self.rootChord+self.tipChord)))+(1/6)*(self.rootChord+self.tipChord-(self.rootChord*self.tipChord)/(self.rootChord+self.tipChord)))

        if self.noseType == "cone":
            Pn = CN_n * (0.6667*self.noseLength)
        elif self.noseType == "ogive":
            Pn = CN_n * (0.466*self.noseLength)
        elif self.noseType == "paraboloid":
            Pn = CN_n * (0.5*self.noseLength)
        elif self.noseType == "ellipsoid":
            Pn = CN_n * (0.3333*self.noseLength)
        else:
            Pn = 0

        netNormalForce = CN_n + CN_f1
        netMoment = Pn + Pf1
        return netMoment/netNormalForce

    

    def temperature(self, altitude):
        if altitude+self.initialAltitude < 11000:
            return self.initialTemperature - self.temperatureLapseRate * (altitude)
        else:
            return self.initialTemperature - self.temperatureLapseRate * (11000-self.initialAltitude)

    def pressure(self, altitude):
        return self.initialPressure * ((self.temperature(altitude)+273.15)/(self.initialTemperature+273.15))**(self.gravity*self.airMolarMass/(self.universeGasConstant*self.temperatureLapseRate))

    def dynamicPressure(self, windSpeed, rocketVelocityX, rocketVelocityY, altitude):
        return 0.5 * self.airDensity(altitude) * self.airSpeed(windSpeed, rocketVelocityX, rocketVelocityY)**2

    def airDensity(self, altitude):
        return self.pressure(altitude)/(self.SpecificGasConstant*(self.temperature(altitude)+273.15))

    def speedOfSound(self, altitude):
        return np.sqrt(self.AdiabaticIndex*self.SpecificGasConstant*(self.temperature(altitude)+273.15))

    def machNumber(self, windSpeed, rocketVelocityX, rocketVelocityY, altitude):
        return self.airSpeed(windSpeed, rocketVelocityX, rocketVelocityY)/self.speedOfSound(altitude)

    def airSpeed(self, windspeed, rocketVelocityX, rocketVelocityY):
        return np.sqrt((rocketVelocityX-windspeed*np.sin(self.windHeading))**2 + (rocketVelocityY**2) + (windspeed*np.cos(self.windHeading))**2)

    def rocketSpeed(self, rocketVelocityX, rocketVelocityY):
        return np.sqrt(rocketVelocityX**2 + rocketVelocityY**2)

    def angle_of_attack(self, theta, x_velocity, y_velocity):
        return np.arctan2(x_velocity, y_velocity) - theta

    def dynamicViscosity(self, altitude):
        return 1.458*10**(-6)*((self.temperature(altitude)+273.15)**(3/2))/(self.temperature(altitude)+273.15+110.4)

    def kinematicViscosity(self, altitude):
        return self.dynamicViscosity(altitude)/self.airDensity(altitude)

    def reynoldsNumber(self, windSpeed, rocketVelocityX, rocketVelocityY, altitude, characteristicLength):
        return self.airSpeed(windSpeed, rocketVelocityX, rocketVelocityY) * characteristicLength / self.kinematicViscosity(altitude)

    def zeroAngleDragCoefficient(self, time, windSpeed, rocketVelocityX, rocketVelocityY, altitude, L):
        return self.C_D_friction(windSpeed, rocketVelocityX, rocketVelocityY, altitude, L) + self.nosePressureDrag(windSpeed, rocketVelocityX, rocketVelocityY, altitude) + self.finPressureDrag(time, windSpeed, rocketVelocityX, rocketVelocityY, altitude) + self.baseDragCoefficient(time, windSpeed, rocketVelocityX, rocketVelocityY, altitude) + self.parasiticDrag(windSpeed, rocketVelocityX, rocketVelocityY, altitude)

    def C_D_friction(self, windSpeed, rocketVelocityX, rocketVelocityY, altitude, L):
        R = self.reynoldsNumber(windSpeed, rocketVelocityX, rocketVelocityY, altitude, L)
        M = self.machNumber(windSpeed, rocketVelocityX, rocketVelocityY, altitude)
        R_crit = 51*(self.surfaceRoughness/L)**-1.039
        if R < 10**4:
            C_f = 1.48*10**-2
        elif R < R_crit:
            C_f = 1/(1.5*np.log(R)-5.6)**2
        else:
            C_f = 0.032*(self.surfaceRoughness/L)**0.2

        C_Mach = (1-0.1*M**2)
        C_f_rough = 0.032*(self.surfaceRoughness/L)**0.2 * C_Mach
        C_f_component = max(C_f, C_f_rough)
        K_body = 1 + 1/(2*self.finenessRatio)
        K_fin = 1 + 2*self.finThickness/self.meanAerodynamicChordLengthOfFin
        C_D_body = C_f_component * K_body * (self.A_wet_body/self.A_ref)
        C_D_fins = C_f_component * K_fin * (self.A_wet_fins/self.A_ref)
        return C_D_body + C_D_fins

    def nosePressureDrag(self, windSpeed, rocketVelocityX, rocketVelocityY, altitude):
        kappa=(1/1) #rho_t/rho in current case both are equal so for performance, it will be simplified.
        M = self.machNumber(windSpeed, rocketVelocityX, rocketVelocityY, altitude)
        epsilon = np.arctan(self.bodyDiameter/(2*self.noseLength))
        gamma = 1.4
        C_D_At_M1 = np.sin(epsilon)
        slope = 4/(gamma+1) * (1-0.5*C_D_At_M1)
        if M < 0.8:
            coneDragCoefficient = 0.8*np.sin(self.jointAngle)**2
        elif M < 1.2:
            coneDragCoefficient = (3*slope+C_D_At_M1-2*np.sin(self.jointAngle)**2)*(M-0.8)+0.8*np.sin(self.jointAngle)**2
        else:
            coneDragCoefficient = slope*M + C_D_At_M1
        return  (0.72*(kappa-0.5)**2 + 0.82)* coneDragCoefficient #the correction factor is only used for ogival shapes, which we are, therefor for now well assume it  

    def finPressureDrag(self, time, windSpeed, rocketVelocityX, rocketVelocityY, altitude):
        M = self.machNumber(windSpeed, rocketVelocityX, rocketVelocityY, altitude)
        if M < 0.9: #rounded edge
            perpindicularLeadingEdgeDrag = (1-M**2)**(-0.417) -1
        elif M < 1:
            perpindicularLeadingEdgeDrag = 1 - 1.785*(M-0.9)
        else:
            perpindicularLeadingEdgeDrag = 1.214 - (0.502/M**2) + (0.1095/M**4)
        leadingEdgeDrag = perpindicularLeadingEdgeDrag * np.cos(self.leadingEdgeAngle)**2 #angled fin, sweep angle
        trailingEdgeDrag = (1/2) * self.baseDragCoefficient(time, windSpeed, rocketVelocityX, rocketVelocityY, altitude) #rounded edge
        C_D_fin = leadingEdgeDrag + trailingEdgeDrag
        A_FIN = self.numberOfFins * self.finThickness * self.finSpan
        return (A_FIN/self.A_ref) * C_D_fin

    def baseDragCoefficient(self, time, windSpeed, rocketVelocityX, rocketVelocityY, altitude):
        if time < self.timeCurve[-1]:
            A_motor = (np.pi/4) * self.motorDiameter ** 2
        else:
            A_motor = 0
        M = self.machNumber(windSpeed, rocketVelocityX, rocketVelocityY, altitude)
        if M < 1:
            CD_base = 0.12 + 0.13*M**2
        else:
            CD_base = 0.25/M
        return CD_base * (self.A_ref-A_motor)/self.A_ref

    def parasiticDrag(self, windSpeed, rocketVelocityX, rocketVelocityY, altitude):
        M = self.machNumber(windSpeed, rocketVelocityX, rocketVelocityY, altitude)
        if M < 1:
            qRatio = 1 + (M**2/4) + (M**4/40)
            C_D_base = 0.12 + 0.13*M**2
        else:
            qRatio = 1.84 - (0.76/M**2) + (0.166/M**4) + (0.035/M**6)
            C_D_base = 0.25/M
        C_D_stag = 0.85 * qRatio
        C_D_parasitic = max(1.3-0.3*(self.launchGuideLength/self.launchGuideOuterDiameter), 1) * C_D_stag
        outer_Area = np.pi * (self.launchGuideOuterDiameter/2)**2
        inner_Area = np.pi * (self.launchGuideInnerDiameter/2)**2
        A_parasitic = outer_Area - inner_Area
        f = 0.25/((np.log10((self.launchGuideRoughness/(3.7*self.launchGuideInnerDiameter)) + (5.74/(self.reynoldsNumber(windSpeed, rocketVelocityX, rocketVelocityY, altitude, self.launchGuideInnerDiameter)**0.9))))**2)
        C_D_tube = f * (self.launchGuideLength/self.launchGuideInnerDiameter)
        return (C_D_tube * inner_Area + 0.7*(C_D_parasitic + C_D_base)*A_parasitic)/self.A_ref

    def axialDragCoefficient(self, AOA):
        if abs(AOA) <= 17:
            scalefunction = (-3/24565)*abs(AOA)**3 + (9/2890)*AOA**2 + 1
        else:
            scalefunction = (13/1945085)*(AOA-90)**3 + (2847/3890170)*(AOA-90)**2
        return self.zeroAngleDragCoefficient() * scalefunction

    def axialDragForce(self):
        rho = self.airDensity()
        V = self.airSpeed()
        C_A = self.axialDragCoefficient
        return (1/2)*rho*V**2 * self.A_ref * C_A

    def normalForceFinCoefficientDerivative(self, AOA, windSpeed, rocketVelocityX, rocketVelocityY, altitude):
        M = self.machNumber(windSpeed, rocketVelocityX, rocketVelocityY, altitude)
        midchordSweepAngle = np.arctan(self.finSpan/(self.sweepLength+self.tipChord/2-self.rootChord/2))
        beta = np.sqrt(abs(M**2-1))
        C_Nalpha0 = 2*np.pi/beta
        AR = 2*(self.finSpan**2)/self.A_fin
        F_D = AR/((1/(2*np.pi))*C_Nalpha0*np.cos(midchordSweepAngle))
        C_N_1fin = (C_Nalpha0 * F_D * (self.A_fin/self.A_ref) * np.cos(midchordSweepAngle))/(2+F_D*np.sqrt(1+4/(F_D**2)))    #(2*np.pi*((finSpan**2)/A_ref))/(2+((beta*finSpan**2)/(A_fin*np.cos(midchordSweepAngle))))
        C_N_fins = (self.numberOfFins/2) * C_N_1fin
        K_TB = 1 + (self.bodyRadius/(self.finSpan+self.bodyRadius))
        return K_TB * C_N_fins

    def windSpeed(self, altitude):
        return self.initWindSpeed * (altitude/self.initialAltitude)


    def normalForceCoefficientDerivative(self, AOA, windSpeed, rocketVelocityX, rocketVelocityY, altitude):
        C_N_finsWithInterference = self.normalForceFinCoefficientDerivative(AOA, windSpeed, rocketVelocityX, rocketVelocityY, altitude)
        C_N_nose = 2*np.cos(AOA)
        return C_N_nose + C_N_finsWithInterference

    def pitchMomentCoefficientDerivative(self, AOA, windSpeed, rocketVelocityX, rocketVelocityY, altitude):
        C_m_nose = (2*np.sin(AOA))/(self.A_ref*self.bodyDiameter)*(self.noseLength*self.A_ref-self.noseVolume)
        M = self.machNumber(windSpeed, rocketVelocityX, rocketVelocityY, altitude)
        AR = 2*self.finSpan**2/self.A_fin
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
            self.meanAerodynamicChordLengthOfFin*((AR*np.sqrt(abs(2**2-1)) - 0.67)/(2*AR*np.sqrt(abs(2**2-1))-1)),   #p(2) = f(2)
            ((0.68*AR)/(np.sqrt(3)*(2*AR*np.sqrt(3)-1)**2)),  #p'(2) = f'(2)
            0,  #p''(2)
            0   #p'''(2)
        ])

        coefficients = np.linalg.solve(A,b)

        if M <= 0.5:
            x_f = (self.sweepLength/3)*((self.rootChord+2*self.tipChord)/(self.rootChord+self.tipChord))+(1/6)*((self.rootChord**2+self.tipChord**2+self.rootChord*self.tipChord)/(self.rootChord+self.tipChord))
        elif M <=2:
            x_f = coefficients[0] + coefficients[1]*M + coefficients[2]*M**2 + coefficients[3]*M**3 + coefficients[4]*M**4 + coefficients[5]*M**5
        else:
            x_f = self.meanAerodynamicChordLengthOfFin*((AR*beta - 0.67)/(2*AR*beta-1))
        x_fin = self.leadingEdgeToNosecone + x_f
        C_m_fins = self.normalForceFinCoefficientDerivative(AOA, windSpeed, rocketVelocityX, rocketVelocityY, altitude) * x_fin/self.bodyDiameter
        return C_m_nose+C_m_fins

    def normalForce(self, AOA, windSpeed, rocketVelocityX, rocketVelocityY, altitude):
        rho = self.airDensity(altitude)
        V = self.airSpeed(windSpeed, rocketVelocityX, rocketVelocityY)
        C_N_alpha = self.normalForceCoefficientDerivative(AOA,windSpeed, rocketVelocityX, rocketVelocityY, altitude)
        return (1/2)*rho*V**2 * self.A_ref * self.bodyDiameter * C_N_alpha * AOA

    def pitchMoment(self, AOA,windSpeed, rocketVelocityX, rocketVelocityY, altitude):
        rho = self.airDensity(altitude)
        V = self.airSpeed(windSpeed,rocketVelocityX,rocketVelocityY)
        C_m_alpha = self.pitchMomentCoefficientDerivative(AOA, windSpeed, rocketVelocityX, rocketVelocityY, altitude)
        return (1/2)*rho*V**2 * self.A_ref * C_m_alpha * AOA

    def centerOfPressure(self, AOA,windSpeed, rocketVelocityX, rocketVelocityY, altitude):
        if AOA != 0:
            return (self.pitchMomentCoefficientDerivative(AOA,windSpeed, rocketVelocityX, rocketVelocityY, altitude)/self.normalForceCoefficientDerivative(AOA,windSpeed, rocketVelocityX, rocketVelocityY, altitude))*self.bodyDiameter
        else:
            return self.CP

    def windSpeedatAltitude(self, altitude):
        return self.initWindSpeed + 0*altitude

    def moment(self, time, theta, rocketVelocityX, rocketVelocityY, altitude):
        AOA = self.angle_of_attack(theta, rocketVelocityX, rocketVelocityY)
        windSpeed = self.windSpeedatAltitude(altitude)
        CG = self.centerOfGravity(time)
        return self.normalForce(AOA, windSpeed, rocketVelocityX, rocketVelocityY, altitude)*(self.centerOfPressure(AOA, windSpeed, rocketVelocityX, rocketVelocityY, altitude)-CG)

    # #Functions






    def moment_of_inertia(self, current_mass, time):
        I_center = (1/12) * current_mass * self.rocketLength ** 2 #assuming uniform density for now
        d = self.centerOfGravity(time) - self.rocketLength/2 #distance between center and CG
        return I_center + current_mass * d**2 #parallel axis theorem

    def normal_force(self, theta, y_pos, x_velocity, y_velocity):
        dynamicPressure = 0.5*self.air_density(y_pos)*self.netSpeed(x_velocity, y_velocity)**2 #Assuming incompressible air for now
        Area = np.pi * self.radius**2
        return dynamicPressure * Area * self.CN_alpha * self.angle_of_attack(theta, x_velocity,y_velocity) 

    def yaw_torque(self, theta, y_pos, x_velocity, y_velocity, time):
        return self.normal_force(theta, y_pos, x_velocity, y_velocity) * (self.CP-self.centerOfGravity(time))

    def alpha(self, time, theta, y_pos, x_velocity, y_velocity):
        return self.moment(time, theta, x_velocity, y_velocity, y_pos) / self.moment_of_inertia(self.mass(time),time)

    def air_density(self, y_pos):
        return 1.225 * (1-((0.0065/288.15)*y_pos))**((self.gravity/(287.05*0.0065))-1)

    def netSpeed(self, x_velocity, y_velocity):
        return np.sqrt((x_velocity**2) + (y_velocity**2))

    def on_rod(self, x_pos, y_pos):
        global leftRod
        if self.leftRod:
            return False
        if np.sqrt((x_pos**2) + (y_pos**2)) <= self.rodLength:
            return True
        leftRod = True
        return False

    def thrust_magnitude(self, time):
        if time > self.timeCurve[-1]:
            return 0
        if time <= self.timeCurve[0]:
            return (self.thrustCurve[1]*time/self.timeCurve[1])
        for i in range(len(self.timeCurve)-1):
            if self.timeCurve[i] < time <= self.timeCurve[i+1]:
                Time1= self.timeCurve[i]
                Time2= self.timeCurve[i+1]
                Thrust1= self.thrustCurve[i]
                Thrust2= self.thrustCurve[i+1]
                return (Thrust1 + (Thrust2 - Thrust1) * (time - Time1) / (Time2 - Time1))
        return 0

    def x_component(self, magnitude, theta):
        return magnitude * np.sin(theta)

    def y_component(self, magnitude, theta):
        return magnitude * np.cos(theta)

    def acceleration_x(self, time, y_pos, x_velocity, y_velocity, theta):
        return self.netForce_x(time, y_pos, x_velocity, y_velocity, theta) / self.mass(time)

    def acceleration_y(self, time, y_pos, x_velocity, y_velocity, theta):
        return self.netForce_y(time, y_pos, x_velocity, y_velocity, theta) / self.mass(time)

    def netForce_x(self, time, y_pos, x_velocity, y_velocity, theta):
        if y_pos == 0:
            return self.x_component(self.thrust_magnitude(time), theta)
        else:
            return self.x_component(self.thrust_magnitude(time), theta) - self.x_component(self.drag_magnitude(y_pos, x_velocity, y_velocity), np.arctan2(x_velocity, y_velocity))    

    def netForce_y(self, time, y_pos, x_velocity, y_velocity, theta):
        if y_pos == 0:
            return self.y_component(self.thrust_magnitude(time), theta)
        else:
            return self.y_component(self.thrust_magnitude(time), theta) - self.weight(time) - self.y_component(self.drag_magnitude(y_pos, x_velocity, y_velocity),np.arctan2(x_velocity, y_velocity))

    def drag_magnitude(self, y_pos, x_velocity, y_velocity):
        speed = self.netSpeed(x_velocity, y_velocity)
        if y_velocity >= 0:
            A, cd = 3.14 * self.radius ** 2, self.dragCoefficient        
        else:
            A, cd = self.parachuteArea, self.parachuteDragCoefficient
        return (0.5 * self.air_density(y_pos) * cd * A * speed**2)

    def weight(self, time):
        return self.mass(time) * self.gravity

    def mass(self, time):
        currentImpulse=0
        netImpulse=0
        for i in range(1,len(self.timeCurve)):
            if time < self.timeCurve[0]:
                currentImpulse = self.thrust_magnitude(time) * time
            elif self.timeCurve[i] <= time:
                currentImpulse += 0.5*(self.thrustCurve[i]+self.thrustCurve[i-1])*(self.timeCurve[i]-self.timeCurve[i-1])
            else:
                currentImpulse += 0.5*(self.thrust_magnitude(time)+self.thrustCurve[i-1]) * (time-self.timeCurve[i-1])
                break
        for i in range(1,len(self.timeCurve)):
            netImpulse += (0.5*(self.thrustCurve[i]+self.thrustCurve[i-1])*(self.timeCurve[i]-self.timeCurve[i-1]))
        currentFuelMass=self.fuelMass*(1-(currentImpulse/netImpulse))
        return self.dryMass + currentFuelMass

    def centerOfGravity(self, time):
        return (self.CG_dry*self.dryMass+self.CT*(self.mass(time)-self.dryMass))/(self.mass(time))
        

    def rk4_step(self, time, x_pos, y_pos, x_velocity, y_velocity, theta, omega):
        k1_x = x_velocity
        k1_y = y_velocity
        k1_vx = self.acceleration_x(time, y_pos, x_velocity, y_velocity, theta)
        k1_vy = self.acceleration_y(time, y_pos, x_velocity, y_velocity, theta)
        k1_0 = omega
        k1_w = self.alpha(time, theta, y_pos, x_velocity, y_velocity)

        k2_x = x_velocity + k1_vx * self.timeStep/2
        k2_y = y_velocity + k1_vy * self.timeStep/2
        k2_vx = self.acceleration_x(time + self.timeStep/2, y_pos + k1_y * self.timeStep/2, x_velocity + k1_vx * self.timeStep/2, y_velocity + k1_vy * self.timeStep/2, theta + k1_0 * self.timeStep/2)
        k2_vy = self.acceleration_y(time + self.timeStep/2, y_pos + k1_y * self.timeStep/2, x_velocity + k1_vx * self.timeStep/2, y_velocity + k1_vy * self.timeStep/2, theta + k1_0 * self.timeStep/2)
        k2_0 = omega + k1_w * self.timeStep/2
        k2_w = self.alpha(time + self.timeStep/2, theta + k1_0 * self.timeStep/2, y_pos + k1_y * self.timeStep/2, x_velocity + k1_vx * self.timeStep/2, y_velocity + k1_vy * self.timeStep/2)

        k3_x = x_velocity + k2_vx * self.timeStep/2
        k3_y = y_velocity + k2_vy * self.timeStep/2
        k3_vx = self.acceleration_x(time + self.timeStep/2, y_pos + k2_y * self.timeStep/2, x_velocity + k2_vx * self.timeStep/2, y_velocity + k2_vy * self.timeStep/2, theta + k2_0 * self.timeStep/2)
        k3_vy = self.acceleration_y(time + self.timeStep/2, y_pos + k2_y * self.timeStep/2, x_velocity + k2_vx * self.timeStep/2, y_velocity + k2_vy * self.timeStep/2, theta + k2_0 * self.timeStep/2)
        k3_0 = omega + k2_w * self.timeStep/2
        k3_w = self.alpha(time + self.timeStep/2, theta + k2_0 * self.timeStep/2, y_pos + k2_y * self.timeStep/2, x_velocity + k2_vx * self.timeStep/2, y_velocity + k2_vy * self.timeStep/2)

        k4_x = x_velocity + k3_vx * self.timeStep
        k4_y = y_velocity + k3_vy * self.timeStep
        k4_vx = self.acceleration_x(time + self.timeStep, y_pos + k3_y * self.timeStep, x_velocity + k3_vx * self.timeStep, y_velocity + k3_vy * self.timeStep, theta + k3_0 * self.timeStep)
        k4_vy = self.acceleration_y(time + self.timeStep, y_pos + k3_y * self.timeStep, x_velocity + k3_vx * self.timeStep, y_velocity + k3_vy * self.timeStep, theta + k3_0 * self.timeStep)
        k4_0 = omega + k3_w * self.timeStep
        k4_w = self.alpha(time + self.timeStep, theta + k3_0 * self.timeStep, y_pos + k3_y * self.timeStep, x_velocity + k3_vx * self.timeStep, y_velocity + k3_vy * self.timeStep)


        new_x_pos = x_pos + (self.timeStep/6) * (k1_x + 2*k2_x + 2*k3_x + k4_x)
        new_y_pos = y_pos + (self.timeStep/6) * (k1_y + 2*k2_y + 2*k3_y + k4_y)
        new_x_velocity = x_velocity + (self.timeStep/6) * (k1_vx + 2*k2_vx + 2*k3_vx + k4_vx)
        new_y_velocity = y_velocity + (self.timeStep/6) * (k1_vy + 2*k2_vy + 2*k3_vy + k4_vy)
        new_theta = theta + (self.timeStep/6) * (k1_0 + 2*k2_0 + 2*k3_0 + k4_0)
        new_omega = omega + (self.timeStep/6) * (k1_w + 2*k2_w + 2*k3_w + k4_w)

        if(self.on_rod(new_x_pos, new_y_pos)):
            new_theta = theta
            new_omega = 0
        
        if new_y_velocity < 0:
            new_theta = theta
            new_omega = 0

        new_theta = (new_theta + np.pi) % (2*np.pi) - np.pi

        return new_x_pos, new_y_pos, new_x_velocity, new_y_velocity, new_theta, new_omega
        
    def run_simulation(self):
        global time, x_pos, y_pos, x_velocity, y_velocity, theta, omega
        while self.y_pos >= -0.01:

            current_mass = self.mass(self.time)
            current_thust = self.thrust_magnitude(self.time)
            current_thrust_x = self.x_component(current_thust, self.theta)
            current_thrust_y = self.y_component(current_thust, self.theta)
            current_weight = current_mass * self.gravity
            current_drag = self.drag_magnitude(self.y_pos, self.x_velocity, self.y_velocity)
            current_drag_x = self.x_component(current_drag, np.arctan2(self.x_velocity, self.y_velocity))
            current_drag_y = self.y_component(current_drag, np.arctan2(self.x_velocity, self.y_velocity))
            current_normalForce = self.normal_force(self.theta, self.y_pos, self.x_velocity, self.y_velocity)
            current_inertia = self.moment_of_inertia(current_mass, self.centerOfGravity(self.time))
            
            if self.y_pos == 0:
                current_netForce_x = current_thrust_x
                current_netForce_y = current_thrust_y 
            else:
                current_netForce_x = current_thrust_x - current_drag_x
                current_netForce_y = current_thrust_y - current_weight - current_drag_y

            if self.on_rod(self.x_pos, self.y_pos):
                current_torque = 0
            else:
                current_torque = self.yaw_torque(self.theta,self.y_pos,self.x_velocity,self.y_velocity,self.time)

            current_acceleration_x = current_netForce_x / current_mass
            current_acceleration_y = current_netForce_y / current_mass
            current_alpha = current_torque/current_inertia

            self.timeList.append(self.time)
            self.x_posList.append(self.x_pos)
            self.y_posList.append(self.y_pos)
            self.x_velocityList.append(self.x_velocity)  
            self.y_velocityList.append(self.y_velocity)
            self.thetaList.append(np.degrees(self.theta))
            self.omegaList.append(self.omega)
            self.alphaList.append(current_alpha)
            self.torqueList.append(current_torque)
            self.inertiaList.append(current_inertia)
            self.normalForceList.append(current_normalForce)
            self.massList.append(current_mass)
            self.x_accelerationList.append(current_acceleration_x)
            self.y_accelerationList.append(current_acceleration_y)
            self.x_netForceList.append(current_netForce_x)
            self.y_netForceList.append(current_netForce_y)
            self.weightList.append(current_weight)
            self.x_thrustList.append(current_thrust_x)
            self.y_thrustList.append(current_thrust_y)
            self.x_dragList.append(current_drag_x)
            self.y_dragList.append(current_drag_y)
            self.densityList.append(self.air_density(self.y_pos))
            new_x_pos, new_y_pos, new_x_velocity, new_y_velocity, new_theta, new_omega = self.rk4_step(self.time, self.x_pos, self.y_pos, self.x_velocity, self.y_velocity, self.theta, self.omega)
            self.x_pos = new_x_pos
            self.y_pos = new_y_pos
            self.x_velocity = new_x_velocity
            self.y_velocity = new_y_velocity
            self.theta = new_theta
            self.omega = new_omega
            self.time += self.timeStep


        for i in range(len(self.timeList)):
            print(f"Time: {self.timeList[i]:.3f} s, \
                Mass: {self.massList[i]:.3f} kg, \
                x_pos: {self.x_posList[i]:.3f} m, \
                Altitude: {self.y_posList[i]:.3f} m, \
                x_vel: {self.x_velocityList[i]:.3f} m/s, \
                Velocity: {self.y_velocityList[i]:.3f} m/s, \
                x_acc: {self.x_accelerationList[i]:.3f} m/s^2 \
                acceleration: {self.y_accelerationList[i]:.3f} m/s^2, \
                x_net: {self.x_netForceList[i]:.3f} N \
                Net Force: {self.y_netForceList[i]:.3f} N, \
                Weight: {self.weightList[i]:.3f} N, \
                x_thrust: {self.x_thrustList[i]:.3f} N \
                Thrust: {self.y_thrustList[i]:.3f} N, \
                x_drag: {self.x_dragList[i]:.3f} N, \
                Drag: {self.y_dragList[i]:.3f} N, \
                theta: {self.thetaList[i]:.3f} degrees, \
                omega: {self.omegaList[i]:.3f} rad/s, \
                alpha: {self.alphaList[i]:.3f} rad/s^2, \
                torque: {self.torqueList[i]:.3f} N*m, \
                inertia: {self.inertiaList[i]:.3f} kg*m^2, \
                normalF: {self.normalForceList[i]:.3f} N, \
                Air Density: {self.densityList[i]:.3f} kg/m^3")
            self.rodVelocity = 0
            for i in range(len(self.y_posList)-1):
                dist1=np.sqrt(self.y_posList[i]**2+self.x_posList[i]**2)
                dist2=np.sqrt(self.y_posList[i+1]**2+self.x_posList[i+1]**2)
                if dist1 < self.rodLength <= dist2:
                    self.rodVelocity = np.sqrt((self.y_velocityList[i]**2 + 2 * self.y_accelerationList[i] * (self.rodLength*np.cos(self.launchAngle) - self.y_posList[i]))+(self.x_velocityList[i]**2 + 2 * self.x_accelerationList[i] * (self.rodLength*np.sin(self.launchAngle) - self.x_posList[i])))
                    break
            print(f"Velocity off rod: {self.rodVelocity:.3f} m/s, Apogee: {max(self.y_posList):.3f} m, Max Velocity: {max(self.y_velocityList):.3f} m/s, Max acceleration: {max(self.y_accelerationList):.3f} m/s^2, Time to Apogee: {self.timeList[self.y_posList.index(max(self.y_posList))]:.3f} s, Flight Time: {self.timeList[-1]:.3f} s, Ground hit velocity: {self.y_velocityList[-1]:.3f} m/s, Range from launch: {self.x_posList[-1]:.3f} m")

    def get_results(self):
        """Return all recorded time-series data as a dict of lists (used to
        build a pandas DataFrame in main.py - no shared file needed)."""
        return {
            "time": self.timeList,
            "mass": self.massList,
            "x_pos": self.x_posList,
            "y_pos": self.y_posList,
            "x_vel": self.x_velocityList,
            "velocity": self.y_velocityList,
            "x_acc": self.x_accelerationList,
            "acceleration": self.y_accelerationList,
            "x_net": self.x_netForceList,
            "net force": self.y_netForceList,
            "weight": self.weightList,
            "x_thrust": self.x_thrustList,
            "thrust": self.y_thrustList,
            "x_drag": self.x_dragList,
            "drag": self.y_dragList,
            "theta": self.thetaList,
            "omega": self.omegaList,
            "alpha": self.alphaList,
            "torque": self.torqueList,
            "inertia": self.inertiaList,
            "normalF": self.normalForceList,
            "air density": self.densityList,
        }

    def results_csv_bytes(self):
        """Build a CSV of the last run's results entirely in memory (no shared
        output.csv on disk, so one session's download never races another's)."""
        results = self.get_results()
        buf = io.StringIO()
        writer = csv.writer(buf)
        writer.writerow(results.keys())
        writer.writerows(zip(*results.values()))
        return buf.getvalue().encode("utf-8")

    def get_summary(self):
        if not self.timeList:
            raise RuntimeError(
                "run_simulation() must be called before get_summary()")

        rodVelocity = 0
        for i in range(len(self.y_posList) - 1):
            dist1 = np.sqrt(self.y_posList[i]**2 + self.x_posList[i]**2)
            dist2 = np.sqrt(self.y_posList[i + 1]**2 +
                            self.x_posList[i + 1]**2)
            if dist1 < self.rodLength <= dist2:
                rodVelocity = np.sqrt(
                    (self.y_velocityList[i]**2 +
                     2 * self.y_accelerationList[i] *
                     (self.rodLength * np.cos(self.launchAngle) -
                      self.y_posList[i])) +
                    (self.x_velocityList[i]**2 +
                     2 * self.x_accelerationList[i] *
                     (self.rodLength * np.sin(self.launchAngle) -
                      self.x_posList[i])))
                break

        apogee = max(self.y_posList)
        return {
            "rod_velocity": rodVelocity,
            "apogee": apogee,
            "max_velocity": max(self.y_velocityList),
            "max_acceleration": max(self.y_accelerationList),
            "time_to_apogee": self.timeList[self.y_posList.index(apogee)],
            "flight_time": self.timeList[-1],
            "ground_hit_velocity": self.y_velocityList[-1],
            "range": self.x_posList[-1],
        }
            

    def plot_results(self):
        fig, ax1 = plt.subplots(figsize=(10, 8))
        ax1.plot(self.timeList, self.y_posList, color='blue', label='Altitude (m)')
        ax1.set_ylabel('Altitude (m)', color='blue')
        ax1.tick_params(axis='y', labelcolor='blue')

        ax2 = ax1.twinx()

        ax2.plot(self.timeList, self.y_velocityList, color='red', label='Velocity (m/s)')
        ax2.set_ylabel('Velocity (m/s)', color='red')
        ax2.tick_params(axis='y', labelcolor='red')
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

        plt.title('Vertical Motion vs Time')
        plt.tight_layout()
        plt.show()

    

test = RocketSimNew()
test.run_simulation()

sys.stdout.flush()
#plot_results()