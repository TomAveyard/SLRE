from PropTools.Utils.constants import G
from PropTools.Utils.mathsUtils import radiusOfCurvature3Points2D
from PropTools.Thermo.fluidDynamics import smoothFrictionFactor
from math import pi, tanh, sqrt

# Gas Side Correlations

def bartzEquation(throatDiameter: float, viscosity: float, specificHeat: float, prandtlNumber: float, chamberPressure: float, chamberTemp: float, cStar: float, localArea: float, gasSideWallTemp: float, machNumber: float, gamma: float, throatAverageRadiusOfCurvature: float = None, C1: float = 0.026) -> float:

        firstBracket = C1 / (throatDiameter ** 0.2)
        secondBracket = (((viscosity ** 0.2) * specificHeat) / (prandtlNumber ** 0.6))
        thirdBracket = ((chamberPressure) / (cStar)) ** 0.8

        # If a throat average radius of curvature is given, the correction is applied
        if throatAverageRadiusOfCurvature != None:
            fourthBracket = (throatDiameter / throatAverageRadiusOfCurvature) ** 0.1
        else:
            fourthBracket = 1

        throatArea = pi * (throatDiameter / 2) ** 2

        fifthBracket = (throatArea / localArea) ** 0.9

        sigmaFirstBracket = ((0.5 * (gasSideWallTemp / chamberTemp) * (1 + (((gamma - 1) / 2) * (machNumber ** 2))) + 0.5) ** (-0.68))
        sigmaSecondBracket = (1 + ((gamma - 1) / 2) * (machNumber ** 2)) ** (-0.12)

        sigma = sigmaFirstBracket * sigmaSecondBracket

        return firstBracket * secondBracket * thirdBracket * fourthBracket * fifthBracket * sigma

# Coolant Side Correlations

def dittusBoelterEquation(reynoldsNumber: float, prandtlNumber: float, C1: float = 0.023) -> float:

    return C1 * (reynoldsNumber ** 0.8) * (prandtlNumber ** 0.4)

def siederTateEquation(reynoldsNumber: float, prandtlNumber: float, bulkViscosity: float, surfaceViscosity: float, C1: float = 0.023) -> float:

    return C1 * (reynoldsNumber ** (4/5)) * (prandtlNumber ** (1/3)) * ((bulkViscosity / surfaceViscosity) ** 0.14)

def superCriticalEquation(specificHeat: float, viscosity: float, prandtlNumber: float, massFlowRate: float, hydraulicDiameter: float, coolantBulkTemp: float, coolantSideWallTemp: float, C1: float = 0.023):

    return ((C1 * specificHeat * (viscosity ** 0.2)) / (prandtlNumber ** (2/3))) * (((massFlowRate * G) ** 0.8) / (hydraulicDiameter ** 0.2)) * ((coolantBulkTemp / coolantSideWallTemp) ** 0.55)

def taylorEquation(reynoldsNumber: float, prandtlNumber: float, coolantBulkTemp: float, coolantSurfaceTemp: float, hydraulicDiameter: float, distance: float, C1: float = 0.023) -> float:

    exponent = 0.57 - 1.59 * (hydraulicDiameter / distance)

    return C1 * (reynoldsNumber ** 0.8) * (prandtlNumber ** 0.4) * ((coolantBulkTemp / coolantSurfaceTemp) ** exponent)

def ruanMengEquation(reynoldsNumber: float, prandtlNumber: float, coolantBulkDensity: float, coolantSurfaceDensity: float, hydraulicDiameter: float, distance: float, C1: float = 0.0069) -> float:

    return C1 * (reynoldsNumber ** 0.8) * (prandtlNumber ** 0.66) * ((coolantBulkDensity / coolantSurfaceDensity) ** 0.43) * (1 + 2.4 * (hydraulicDiameter / distance))

def nusseltNumberToHeatTransferCoefficient(nusseltNumber: float, charactersiticLength: float, thermalConductivity: float) -> float:

    return (nusseltNumber * thermalConductivity) / charactersiticLength

# Fin heat transfer

def finEffectiveness(heatTransferCoeffcient: float, ribThickness: float, thermalConductivity: float, channelHeight: float) -> float:
        
    common = sqrt((2 * heatTransferCoeffcient * ribThickness) / thermalConductivity) * (channelHeight / ribThickness)

    return tanh(common) / common

def applyFinEffectiveness(heatTransferCoefficient: float, finCorrectionFactor: float, channelWidth: float, channelHeight: float, ribThickness: float) -> float:

    numerator = channelWidth + (2 * finCorrectionFactor * channelHeight)
    denominator = channelWidth + ribThickness

    return heatTransferCoefficient * (numerator / denominator)

# Curvature correction

def curvatureCorrectionFactor(p1: list, p2: list, p3: list, reynoldsNumber: float, hydraulicDiameter: float) -> float:

    radiusOfCurvature, centrePoint = radiusOfCurvature3Points2D(p1, p2, p3, returnCenterPoint=True)

    if radiusOfCurvature == None:

        return 1

    if centrePoint[1] > p2[1]:
        exponent = 0.05
    else:
        exponent = -0.05

    return (reynoldsNumber * ((hydraulicDiameter / (2 * radiusOfCurvature)) ** 2)) ** exponent

# Roughness correction

def roughnessCorrectionFactor(frictionFactor: float, reynoldsNumber: float, prandtlNumber: float) -> float:

    roughSmoothRatio = frictionFactor / smoothFrictionFactor(reynoldsNumber)

    numerator = 1 + 1.5 * (prandtlNumber ** (-1/6)) * (reynoldsNumber ** (-1/8)) * (prandtlNumber - 1)
    denominator = 1 + 1.5 * (prandtlNumber ** (-1/6)) * (reynoldsNumber ** (-1/8)) * ((prandtlNumber * roughSmoothRatio) - 1)

    return (numerator / denominator) * roughSmoothRatio


