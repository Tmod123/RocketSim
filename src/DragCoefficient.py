import numpy as np
import random

#rocketVariables
referenceViscosity = 1.716*(10**(-5))
referenceTemperature = 273.15
SutherlandConstant = 110.4

cone = "ogive"
L_N = 0.1
L_BT = 0.3
d = 0.025
t=0.002
rootC=0.0508
tipC=0.0508
finLength = 0.03
numberOfFins = 3
S_FF = finLength*rootC*tipC
L_LL = 0.035
LLoutD = 0.007
LLinD = 0.005

#functions
def reynoldsNumber(airDensity, velocity, tempF, characteristicLength): 
    return airDensity * velocity * characteristicLength / dynamicViscosity(tempF)

def dynamicViscosity(tempF):
    tempK = (tempF-32)*(5/9)+273.15
    return referenceViscosity * ((tempK/referenceTemperature)**(3/2)) * ((referenceTemperature+SutherlandConstant)/(tempK + SutherlandConstant))

def skinFrictionCoefficient(airDensity, velocity, tempF, characteristicLength):
    Re = reynoldsNumber(airDensity, velocity, tempF, characteristicLength)
    if Re <= 10**5 or random.randint(0,1) == 0:
        return 0.664/np.sqrt(Re)
    else:
        return 0.0576/(Re**0.2)
    
def wettedSurfaceToCrossSectionRatio():
    if cone == "conical":
        ratio_N = 2*(L_N/d)*np.sqrt(1+(1/(4*((L_N/d)**2))))
    elif cone == "ogive":
        ratio_N = (4/np.pi)*((L_N/d)**2)*np.arcsin(1/(2*L_N/d))
    
    ratio_BT = 4*L_BT/d
    return ratio_N + ratio_BT

def partialSum(airDensity, velocity, tempF):
    C_f = skinFrictionCoefficient(airDensity, velocity, tempF, L_N+L_BT)
    S_W2S_BT = wettedSurfaceToCrossSectionRatio()
    return 1.02*C_f*(1+1.5/(((L_BT+L_N)/d)**(3/2)))*(S_W2S_BT)


def baseDragCoefficient(partialCDSum):
    return 0.029/np.sqrt(partialCDSum)

def finDragSurfaceArea(airDensity, velocity, tempF):
    C_f = skinFrictionCoefficient(airDensity, velocity, tempF, rootC)
    return 2*C_f*(1+2*(t/rootC))


def finDragCoefficientAtZero(airDensity, velocity, tempF):
    fDSA = finDragSurfaceArea(airDensity, velocity, tempF)
    S_F = numberOfFins * S_FF
    S_BT = (np.pi/4)*(d)**2
    return fDSA * S_F/S_BT

def finAndBodyInterferenceCoefficient(airDensity, velocity, tempF):
    fDSA = finDragSurfaceArea(airDensity, velocity, tempF)
    S_BT = (np.pi/4)*(d)**2
    return fDSA * (rootC/S_BT) * (d/2) * numberOfFins

def launchLugDragCoefficient():
    S_LL_W = np.pi*L_LL*(LLoutD+LLinD)
    S_LL = (np.pi/4)*((LLoutD**2)-(LLinD)**2)
    S_BT = (np.pi/4)*(d)**2
    return (1.2*S_LL + 0.045*S_LL_W)/S_BT

def netDragCoefficientAtZero(airDensity, velocity, tempF):
    partialCDSum = partialSum(airDensity, velocity, tempF)
    C_D_B = baseDragCoefficient(partialCDSum)
    C_D_0F = finDragCoefficientAtZero(airDensity, velocity, tempF)
    C_D_int = finAndBodyInterferenceCoefficient(airDensity, velocity, tempF)
    C_D_LL = launchLugDragCoefficient()

    return partialCDSum + C_D_B + C_D_0F + C_D_int + C_D_LL


if __name__ == "__main__":
    #print(reynoldsNumber(1.225, 100, 59, 0.4))
    print(netDragCoefficientAtZero(1.225, 100, 59))