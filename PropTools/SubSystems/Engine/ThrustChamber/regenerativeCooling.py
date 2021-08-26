from matplotlib.pyplot import cool
from PropTools.SubSystems.Engine.Propellant.propellant import Propellant
import numpy as np
from math import sqrt, pi
from PropTools.SubSystems.Engine.ThrustChamber.thrustChamber import ThrustChamber
from PropTools.Utils.constants import G
from PropTools.Utils.mathsUtils import radiusToArea, distanceBetweenTwoPoints
import matplotlib.pyplot as plt

# Class to store information about the cooling channel design
# Channels are modelled as a sector of a annulus
class CoolingChannels:

    def __init__(self,
        massFlowRate, 
        numberOfChannels,
        wallThickness,
        midRibThickness,
        channelHeight, 
        wallConductivity, 
        wallRoughnessHeight=0):

        self.massFlowRate = massFlowRate
        self.numberOfChannels = numberOfChannels
        self.wallThickness = wallThickness
        self.midRibThickness = midRibThickness
        self.channelHeight = channelHeight
        self.wallConductivity = wallConductivity
        self.wallRoughnessHeight = wallRoughnessHeight

        self.channelInstance = ChannelDimensions(self)


# Class to store the dimensions of a single instance of the cooling channels at a certain area
class ChannelDimensions:

    def __init__(self, coolingChannels: CoolingChannels):

        self.coolingChannels = coolingChannels

        self.bottomWidth = None
        self.midWidth = None
        self.topWidth = None
        self.bottomRadius = None
        self.midRadius = None
        self.topRadius = None
        self.channelAngle = None
        self.ribAngle = None
        self.thrustChamberRadius = None
        self.individualChannelArea = None
        self.totalChannelArea = None

    def getChannelDimensions(self, thrustChamberArea):
        
        self.thrustChamberRadius = sqrt(thrustChamberArea / pi)

        self.bottomRadius = self.thrustChamberRadius + self.coolingChannels.wallThickness
        self.topRadius = self.bottomRadius + self.coolingChannels.channelHeight
        self.midRadius = (self.bottomRadius + self.topRadius) / 2
        self.height = self.topRadius - self.bottomRadius

        self.ribAngle = (360 * self.coolingChannels.midRibThickness) / (2 * pi * self.midRadius)
        self.channelAngle = (360 - self.ribAngle * self.coolingChannels.numberOfChannels) / self.coolingChannels.numberOfChannels

        self.bottomWidth = 2 * pi * self.bottomRadius * self.channelAngle / 360
        self.topWidth = 2 * pi * self.topRadius * self.channelAngle / 360
        self.midWidth = 2 * pi * self.midRadius * self.channelAngle / 360

        self.individualChannelArea = 0.5 * self.height * (self.bottomWidth + self.topWidth)
        self.totalChannelArea = self.individualChannelArea * self.coolingChannels.numberOfChannels
   
# Class to store information about the regenerative cooling system
# Performs an iterative calculation for every point defined by the input thrust chamber object
# with the channel dimensions defined by the cooling channel object
class RegenerativeCooling:

    def __init__(self, thrustChamber: ThrustChamber, coolingChannels: CoolingChannels, coolantInletState: Propellant):

        self.thrustChamber = thrustChamber
        self.coolingChannels = coolingChannels
        self.coolantInletState = coolantInletState

    def calculate(self, convergenceCriteria=0.01):

        i = len(self.thrustChamber.axialCoords)

        # Initialise arrays to store relevant info about the heat transfer
        self.heatFluxes = np.zeros(i)
        self.adiabaticWallTemps = np.zeros(i)
        self.gasSideWallTemps = np.zeros(i)
        self.coolantSideWallTemps = np.zeros(i)
        self.coolantBulkTemps = np.zeros(i)
        self.coolantPressures = np.zeros(i)

        # -1 from the length to get the last index of the array
        i -= 1

        # Get values for first station at very end of nozzle
        self.coolantBulkTemps[i] = self.coolantInletState.T
        self.coolantPressures[i] = self.coolantInletState.P
        self.adiabaticWallTemps[i] = self.getLocalAdiabaticWallTempNozzle(self.thrustChamber.exitArea)

        # Propellant objects for calculating the changes in the coolant state as it passes through the channels
        coolantState = Propellant(self.coolantInletState.name)
        coolantState.defineState("T", self.coolantInletState.T, "P", self.coolantInletState.P)
        surfaceCoolantState = Propellant(self.coolantInletState.name)

        i -= 1

        # Initialise variables for loop
        area = 0
        adiabaticWallTemp = 0
        localMachNumber = 0
        gasSideHeatTransferCoefficient = 0
        gasSideHeatFlux = 0
        newGasSideWallTemp = convergenceCriteria * 10 # So the convergence criteria isn't met before the loop begins
        coolantSideWallTemp = 0
        coolantSideWallTempPrev = coolantState.T
        coolantVelocity = 0
        coolantReynoldsNumber = 0
        coolantPrandtlNumber = 0
        coolantNusseltNumber = 0
        coolantSideHeatTransferCoefficient = 0
        
        totalIterations = 0

        while i > 0:

            iterations = 0

            # Get the area, adiabatic wall temperature, and gas mach number for the station before the loop,
            # as these are constant for a station
            area = radiusToArea(self.thrustChamber.radialCoords[i])

            if i > 0:
                localMachNumber = self.getLocalMachNumberNozzle(area)
                adiabaticWallTemp = self.getLocalAdiabaticWallTempNozzle(area)
            elif i < 0:
                localMachNumber = self.getLocalMachNumberChamber(area)
                adiabaticWallTemp = self.getLocalAdiabaticWallTempChamber(area)

            # Guess an initial value for the gas side wall temp
            gasSideWallTemp = self.gasSideWallTemps[i-1]

            while abs(gasSideWallTemp - newGasSideWallTemp) > convergenceCriteria:
                
                # Calculate the gas side heat transfer coefficient and heat flux for the guess value
                gasSideHeatTransferCoefficient = self.bartzEquation(area, gasSideWallTemp, localMachNumber)
                gasSideHeatFlux = gasSideHeatTransferCoefficient * (adiabaticWallTemp - gasSideWallTemp)

                # Calculate flow properties for the channel at the station
                self.coolingChannels.channelInstance.getChannelDimensions(area)
                coolantVelocity = self.massFlowRateToVelocity(self.coolingChannels.massFlowRate, coolantState.D, self.coolingChannels.channelInstance.totalChannelArea)
                coolantReynoldsNumber = self.reynoldsNumber(coolantState.D, coolantVelocity, self.coolingChannels.channelInstance.midWidth, coolantState.viscosity)
                coolantPrandtlNumber = self.prandtlNumber(coolantState.cp, coolantState.viscosity, coolantState.thermalConductivity)

                # Initialise the state of the coolant at the station using the value of the previous station
                coolantSideWallTemp = coolantSideWallTempPrev
                surfaceCoolantState.defineState("T", coolantSideWallTemp, "P", coolantState.P)

                # Reset previous coolant side wall temperature to zero to avoid breaking the loop immediately
                coolantSideWallTempPrev = 0

                # Iterate to find the coolant side wall temperature as the surface viscosity is dependent on the surface temperature
                while abs(coolantSideWallTemp - coolantSideWallTempPrev) > convergenceCriteria:

                    coolantSideWallTempPrev = coolantSideWallTemp

                    # Calculate heat transfer coefficient for the coolant side
                    coolantNusseltNumber = self.siederTateEquation(coolantReynoldsNumber, coolantPrandtlNumber, coolantState.viscosity, surfaceCoolantState.viscosity)
                    coolantSideHeatTransferCoefficient = self.nusseltNumberToHeatTransferCoefficient(coolantNusseltNumber, self.coolingChannels.channelInstance.midWidth, coolantState.thermalConductivity)

                    # Calculate the coolant side wall temperature using the heat flux calculated from the gas side
                    coolantSideWallTemp = (gasSideHeatFlux / coolantSideHeatTransferCoefficient) + coolantState.T
                    surfaceCoolantState.defineState("T", coolantSideWallTemp, "P", coolantState.P)

                # Set coolantSideWallTempPrev to the coverged value so that it can be used to initialise the next station
                # for decreased solve time
                coolantSideWallTempPrev = coolantSideWallTemp

                # Calculate the gas side wall temperature using the gas side heat flux and the coolant side wall temperature
                newGasSideWallTemp = gasSideHeatFlux * (self.coolingChannels.wallThickness / self.coolingChannels.wallConductivity) + coolantSideWallTemp

                # Set the gas side wall temperature guess value to the average of the previous guess value, and the calculated value above
                gasSideWallTemp = (gasSideWallTemp + newGasSideWallTemp) / 2

                iterations += 1

            totalIterations += iterations

            print("Solved station " + str(i) + " in " + str(iterations) + " iterations")
                
            self.heatFluxes[i] = gasSideHeatFlux
            self.adiabaticWallTemps[i] = adiabaticWallTemp
            self.gasSideWallTemps[i] = gasSideWallTemp
            self.coolantSideWallTemps[i] = coolantSideWallTemp
            self.coolantBulkTemps[i] = coolantState.T
            self.coolantPressures[i] = coolantState.P

            stationLength = distanceBetweenTwoPoints([self.thrustChamber.axialCoords[i], self.thrustChamber.radialCoords[i]],
                                                    [self.thrustChamber.axialCoords[i+1], self.thrustChamber.radialCoords[i+1]])
            newCoolantBulkTemp = coolantState.T + ((gasSideHeatFlux * 2 * pi * self.thrustChamber.radialCoords[i] * stationLength) / (self.coolingChannels.massFlowRate * coolantState.cp))
            newCoolantPressure = coolantState.P

            coolantState.defineState("T", newCoolantBulkTemp, "P", newCoolantPressure)

            i -= 1

        print("\nHeat transfer calculaton complete")
        print("Total iterations: " + str(totalIterations))

    def bartzEquation(self, localArea, gasSideWallTemp, localMachNumber):

        firstBracket = (0.026 / self.thrustChamber.throatDiameter)
        secondBracket = (((self.thrustChamber.chamberViscosity ** 0.2) * self.thrustChamber.chamberHeatCapacity) / (self.thrustChamber.chamberPrandtlNumber ** 0.6))
        thirdBracket = ((self.thrustChamber.chamberPressure * G) / self.thrustChamber.cStar) ** 0.8
        fourthBracket = (self.thrustChamber.throatDiameter / self.thrustChamber.throatAverageRadiusOfCurvature) ** 0.1
        fifthBracket = (self.thrustChamber.throatArea / localArea) ** 0.9
        sigma = self.bartzCorrectionFactor(gasSideWallTemp, localMachNumber)

        return firstBracket * secondBracket * thirdBracket * fourthBracket * fifthBracket * sigma

    def bartzCorrectionFactor(self, gasSideWallTemp, localMachNumber):

        denominatorFirstBracket = ((0.5 * (gasSideWallTemp / self.thrustChamber.chamberTemp) * (1 + ((self.thrustChamber.chamberGamma - 1) / 2) * (localMachNumber ** 2)) + 0.5) ** 0.68)
        denominatorSecondBracket = (1 + ((self.thrustChamber.chamberGamma - 1) / 2) * (localMachNumber ** 2)) ** 0.12

        return 1 / (denominatorFirstBracket * denominatorSecondBracket)

    def reynoldsNumber(self, density, velocity, characteristicLength, visocity):

        return (density * velocity * characteristicLength) / visocity

    def prandtlNumber(self, specificHeat, viscosity, thermalConductivity):

        return (specificHeat * viscosity) / thermalConductivity

    def massFlowRateToVelocity(self, massFlowRate, density, area):

        return massFlowRate / (density * area)

    def siederTateEquation(self, reynoldsNumber, prandtlNumber, bulkViscosity, surfaceViscosity, C1=0.023):

        return C1 * (reynoldsNumber ** (4/5)) * (prandtlNumber ** (1/3)) * ((bulkViscosity / surfaceViscosity) ** 0.14)

    def nusseltNumberToHeatTransferCoefficient(self, nusseltNumber, charactersiticLength, thermalConductivity):

        return (nusseltNumber * thermalConductivity) / charactersiticLength

    def getLocalMachNumberNozzle(self, area):

        CEA = self.thrustChamber.getCEAObject()

        expansionRatioAtArea = area / self.thrustChamber.throatArea

        c, ct = self.thrustChamber.convertConditionToCEAFlag(self.thrustChamber.condition, self.thrustChamber.throatCondition)

        return CEA.get_MachNumber(Pc=self.thrustChamber.injectionPressure, MR=self.thrustChamber.mixtureRatio, eps=expansionRatioAtArea, frozen=c, frozenAtThroat=ct)

    # Gets local mach number upstream of the throat as it is not possible to use CEA to get this
    # Uses the Newton-Raphson method to iteratively solve for mach number from area ratio to the specified accuracy
    # Method detailed here: https://www.grc.nasa.gov/www/winddocs/utilities/b4wind_guide/mach.html
    def getLocalMachNumberChamber(self, area, accuracy=0.001):

        gamma = self.thrustChamber.chamberGamma

        R = (area / self.thrustChamber.throatArea) ** 2
        E = (gamma + 1) / (gamma - 1)
        P = 2 / (gamma + 1)
        Q = 1 / E

        a = P ** (1/Q)
        r = (R - 1) / (2 * a)

        X = 1 / (1 + r + sqrt(r * (r + 2)))

        def newtonRaphson(X, P, Q):

            return (P * (X - 1)) / (1 - R * (P + Q * X) ** (-P / Q))

        Xnew = newtonRaphson(X, P, Q)

        while abs(Xnew - X) > accuracy:

            X = Xnew
            Xnew = newtonRaphson(X, P, Q)

        subsonicMachNumber = sqrt(Xnew)

        return subsonicMachNumber

    def getLocalGasTempNozzle(self, area):

        CEA = self.thrustChamber.getCEAObject()

        expansionRatioAtArea = area / self.thrustChamber.throatArea

        c, tc = self.thrustChamber.convertConditionToCEAFlag(self.thrustChamber.condition, self.thrustChamber.throatCondition)

        return CEA.get_Temperatures(Pc=self.thrustChamber.injectionPressure, MR=self.thrustChamber.mixtureRatio, eps=expansionRatioAtArea, frozen=c, frozenAtThroat=tc)[-1]

    def getLocalGasTempChamber(self, area):

        mach = self.getLocalMachNumberChamber(area)
        gamma = self.getLocalGamma(area)

        return self.thrustChamber.chamberTemp / (1 + ((gamma - 1) / 2) * (mach ** 2))

    def getLocalRecoveryFactor(self, area):

        expansionRatioAtArea = area / self.thrustChamber.throatArea

        prandtlNumber = self.thrustChamber.getExitTransportPropertiesAtExpansionRatio(expansionRatioAtArea)[3]

        return prandtlNumber ** 0.33

    def getLocalGamma(self, area):

        CEA = self.thrustChamber.getCEAObject()

        expansionRatioAtArea = area / self.thrustChamber.throatArea

        return CEA.get_exit_MolWt_gamma(Pc=self.thrustChamber.chamberPressure, MR=self.thrustChamber.mixtureRatio, eps=expansionRatioAtArea)[1]

    def getLocalAdiabaticWallTempNozzle(self, area):

        mach = self.getLocalMachNumberNozzle(area)
        r = self.getLocalRecoveryFactor(area)
        gamma = self.getLocalGamma(area)
        gasTemp = self.getLocalGasTempNozzle(area)

        numerator = 1 + (r * ((gamma-1) / 2) * (mach ** 2))
        denominator = 1 + ((gamma-1) / 2) * (mach ** 2)

        return gasTemp * (numerator / denominator)

    def getLocalAdiabaticWallTempChamber(self, area):

        mach = self.getLocalMachNumberChamber(area)
        r = self.thrustChamber.chamberPrandtlNumber ** 0.33
        gamma = self.thrustChamber.chamberGamma
        gasTemp = self.getLocalGasTempChamber(area)

        numerator = 1 + (r * ((gamma-1) / 2) * (mach ** 2))
        denominator = 1 + ((gamma-1) / 2) * (mach ** 2)

        return gasTemp * (numerator / denominator)

testThrustChamber = ThrustChamber('ethanol', 'oxygen', 10*10**3, 25, fac=True, CR=5, ambientPressure=1)

testThrustChamber.getChamberGeometry(1.05,
                                     0.05, 
                                     entranceRadiusOfCurvatureFactor=0.75, 
                                     throatEntranceStartAngle=-135, 
                                     numberOfPointsConverging=50,
                                     numberOfPointsStraight=10)

testThrustChamber.getRaoBellNozzleGeometry(0.6, numberOfPoints=50)
testThrustChamber.getThrustChamberCoords()

testCoolingChannels = CoolingChannels(testThrustChamber.fuelMassFlowRate, 80, 0.001, 0.001, 0.01, 300)

inlet = Propellant(testThrustChamber.fuel.name)
inlet.defineState("T", 298, "P", 40*10**5)

testRegenCooling = RegenerativeCooling(testThrustChamber, testCoolingChannels, inlet)
testRegenCooling.calculate(convergenceCriteria=0.1)

fig, ax = plt.subplots()
ax.plot(testThrustChamber.axialCoords[1:-1], testRegenCooling.heatFluxes[1:-1])
plt.show()

