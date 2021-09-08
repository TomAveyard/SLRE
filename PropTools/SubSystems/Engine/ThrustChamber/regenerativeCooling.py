from thermo import viscosity
from PropTools.SubSystems.Engine.Propellant.propellant import Propellant
import numpy as np
from math import sqrt, pi, log10, tanh
from PropTools.SubSystems.Engine.ThrustChamber.thrustChamber import ThrustChamber
from PropTools.Utils.constants import G
from PropTools.Utils.mathsUtils import radiusToArea, distanceBetweenTwoPoints, radiusOfCurvature3Points2D
from PropTools.SubSystems.Engine.Cycle.component import Component

# Class to store information about the cooling channel design
# Channels are modelled as a sector of a annulus
class CoolingChannels:

    def __init__(self,
        numberOfChannels,
        wallThickness,
        midRibThickness,
        channelHeight, 
        wallConductivity, 
        wallRoughnessHeight):

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
        self.totalChannelAngle = None
        self.individualChannelAngle = None
        self.totalRibAngle = None
        self.individualRibAngle = None
        self.thrustChamberRadius = None
        self.individualChannelArea = None
        self.totalChannelArea = None
        self.sideLength = None
        self.wettedPerimeter = None
        self.hydraulicDiameter = None

    def getChannelDimensions(self, thrustChamberArea):
        
        self.thrustChamberRadius = sqrt(thrustChamberArea / pi)

        self.bottomRadius = self.thrustChamberRadius + self.coolingChannels.wallThickness
        self.topRadius = self.bottomRadius + self.coolingChannels.channelHeight
        self.midRadius = (self.bottomRadius + self.topRadius) / 2
        self.height = self.topRadius - self.bottomRadius

        self.totalRibAngle = (360 * self.coolingChannels.midRibThickness) / (2 * pi * self.midRadius)
        self.individualRibAngle = self.totalRibAngle / self.coolingChannels.numberOfChannels
        self.totalChannelAngle = 360 - self.totalRibAngle
        self.individualChannelAngle = self.totalChannelAngle / self.coolingChannels.numberOfChannels

        self.bottomWidth = 2 * pi * self.bottomRadius * (self.individualChannelAngle / 360)
        self.topWidth = 2 * pi * self.topRadius * (self.individualChannelAngle / 360)
        self.midWidth = 2 * pi * self.midRadius * (self.individualChannelAngle / 360)

        self.individualChannelArea = 0.5 * self.height * (self.bottomWidth + self.topWidth)
        self.totalChannelArea = self.individualChannelArea * self.coolingChannels.numberOfChannels

        self.sideLength = sqrt((((self.topWidth - self.bottomWidth) / 2) ** 2) + (self.height ** 2))
        self.wettedPerimeter = self.bottomWidth + self.topWidth + (self.sideLength * 2)
        self.hydraulicDiameter = 4 * self.individualChannelArea / self.wettedPerimeter

# Class to store information about the regenerative cooling system
# Performs an iterative calculation for every point defined by the input thrust chamber object
# with the channel dimensions defined by the cooling channel object
class RegenerativeCoolingNew(Component):

    def __init__(self, thrustChamber: ThrustChamber, coolingChannels: CoolingChannels, coolantSideHeatTransferCorrelation="Sieder-Tate", includeFinCorrection=True, includeCurvatureCorrection=True, includeRoughnessCorrection=True):

        super().__init__()
        self.type = "regenerative cooling"

        self.thrustChamber = thrustChamber
        self.coolingChannels = coolingChannels
        self.coolantSideHeatTransferCorrelation = coolantSideHeatTransferCorrelation

        self.includeFinCorrection = includeFinCorrection
        self.includeCurvatureCorrection = includeCurvatureCorrection
        self.includeRoughnessCorrection = includeRoughnessCorrection

        self.numberOfStations = len(self.thrustChamber.axialCoords)

        self.heatFluxes = np.zeros(self.numberOfStations)
        self.adiabaticWallTemps = np.zeros(self.numberOfStations)
        self.gasSideWallTemps = np.zeros(self.numberOfStations)
        self.coolantSideWallTemps = np.zeros(self.numberOfStations)
        self.coolantBulkTemps = np.zeros(self.numberOfStations)
        self.coolantPressures = np.zeros(self.numberOfStations)
        self.coolantReynoldsNumbers = np.zeros(self.numberOfStations)
        self.coolantNusseltNumbers = np.zeros(self.numberOfStations)
        self.coolantPrandtlNumbers = np.zeros(self.numberOfStations)
        self.gasSideHeatTransferCoefficients = np.zeros(self.numberOfStations)
        self.coolantSideHeatTransferCoefficients = np.zeros(self.numberOfStations)
        self.channelAreas = np.zeros(self.numberOfStations)

        self.outletState = None
        self.enthalpyChange = None
        self.totalHeatPower = None

        self.inletState = None
        self.massFlowRate = None

    def calculate(self, inletState: Propellant, massFlowRate, convergenceCriteria=0.01):

        print("---")
        print("Starting Heat Transfer Calculation")
        print("Number Of Stations: " + str(self.numberOfStations))






# Class to store information about the regenerative cooling system
# Performs an iterative calculation for every point defined by the input thrust chamber object
# with the channel dimensions defined by the cooling channel object
class RegenerativeCooling(Component):

    def __init__(self, thrustChamber: ThrustChamber, coolingChannels: CoolingChannels, coolantSideHeatTransferCorrelation="Sieder-Tate", includeFinCorrection=True, includeCurvatureCorrection=True, includeRoughnessCorrection=True):

        super().__init__()
        self.type = "regenerative cooling"

        self.thrustChamber = thrustChamber
        self.coolingChannels = coolingChannels
        self.coolantSideHeatTransferCorrelation = coolantSideHeatTransferCorrelation

        self.includeFinCorrection = includeFinCorrection
        self.includeCurvatureCorrection = includeCurvatureCorrection
        self.includeRoughnessCorrection = includeRoughnessCorrection

        self.numberOfStations = len(self.thrustChamber.axialCoords)

        self.heatFluxes = np.zeros(self.numberOfStations)
        self.adiabaticWallTemps = np.zeros(self.numberOfStations)
        self.gasSideWallTemps = np.zeros(self.numberOfStations)
        self.coolantSideWallTemps = np.zeros(self.numberOfStations)
        self.coolantBulkTemps = np.zeros(self.numberOfStations)
        self.coolantPressures = np.zeros(self.numberOfStations)
        self.coolantReynoldsNumbers = np.zeros(self.numberOfStations)
        self.coolantNusseltNumbers = np.zeros(self.numberOfStations)
        self.coolantPrandtlNumbers = np.zeros(self.numberOfStations)
        self.gasSideHeatTransferCoefficients = np.zeros(self.numberOfStations)
        self.coolantSideHeatTransferCoefficients = np.zeros(self.numberOfStations)
        self.channelAreas = np.zeros(self.numberOfStations)

        self.outletState = None
        self.enthalpyChange = None
        self.totalHeatPower = None

        self.inletState = None
        self.massFlowRate = None

    def calculate(self, inletState: Propellant, massFlowRate, convergenceCriteria=0.01):

        print("---")
        print("Starting Heat Transfer Calculation")
        print("Number Of Stations: " + str(self.numberOfStations))

        self.inletState = inletState
        self.outletState = Propellant(self.inletState.name)
        self.massFlowRate = massFlowRate

        i = self.numberOfStations - 1

        # Get values for first station at very end of nozzle
        self.coolantBulkTemps[i] = self.inletState.T
        self.coolantPressures[i] = self.inletState.P
        self.adiabaticWallTemps[i] = self.getLocalAdiabaticWallTempNozzle(self.thrustChamber.exitArea)
        self.gasSideWallTemps[i] = self.adiabaticWallTemps[i] * 0.8

        # Propellant objects for calculating the changes in the coolant state as it passes through the channels
        coolantState = Propellant(self.inletState.name)
        coolantState.defineState("T", self.inletState.T, "P", self.inletState.P)
        surfaceCoolantState = Propellant(self.inletState.name)

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

            if self.thrustChamber.axialCoords[i] >= 0:
                localMachNumber = self.getLocalMachNumberNozzle(area)
                adiabaticWallTemp = self.getLocalAdiabaticWallTempNozzle(area)
            elif self.thrustChamber.axialCoords[i] < 0:
                localMachNumber = self.getLocalMachNumberChamber(area)
                adiabaticWallTemp = self.getLocalAdiabaticWallTempChamber(area)

            # Guess an initial value for the gas side wall temp
            gasSideWallTemp = self.gasSideWallTemps[i+1]

            newGasSideWallTemp = 0

            while abs(gasSideWallTemp - newGasSideWallTemp) > convergenceCriteria:
                
                if self.thrustChamber.axialCoords[i] >= 0:
                    section = "nozzle"
                else:
                    section = "chamber"

                # Calculate the gas side heat transfer coefficient and heat flux for the guess value
                gasSideHeatTransferCoefficient = self.bartzEquation(area, gasSideWallTemp, localMachNumber, section)
                gasSideHeatFlux = gasSideHeatTransferCoefficient * (adiabaticWallTemp - gasSideWallTemp)

                # Calculate flow properties for the channel at the station
                self.coolingChannels.channelInstance.getChannelDimensions(area)
                coolantVelocity = self.massFlowRateToVelocity(self.massFlowRate, coolantState.D, self.coolingChannels.channelInstance.totalChannelArea)
                coolantReynoldsNumber = self.reynoldsNumber(coolantState.D, coolantVelocity, self.coolingChannels.channelInstance.hydraulicDiameter, coolantState.viscosity)
                coolantPrandtlNumber = self.prandtlNumber(coolantState.cp, coolantState.viscosity, coolantState.thermalConductivity)

                # Calculate friction factor
                frictionFactor = self.colebrookEquation(self.coolingChannels.wallRoughnessHeight, self.coolingChannels.channelInstance.hydraulicDiameter, coolantReynoldsNumber, convergenceCriteria=convergenceCriteria)

                if self.includeRoughnessCorrection:
                    smoothFrictionFactor = self.smoothFrictionFactor(coolantReynoldsNumber)
                    roughSmoothRatio = frictionFactor / smoothFrictionFactor

                if self.coolantSideHeatTransferCorrelation == "super critical":

                    # Calculate curvature correction factor if included
                    if self.includeCurvatureCorrection:
                        curvatureCorrectionFactor = self.getCurvatureCorrectionFactor(i, coolantReynoldsNumber, self.coolingChannels.channelInstance.hydraulicDiameter)
                    else:
                        curvatureCorrectionFactor = 1

                    # Calculate roughness correction factor if included
                    if self.includeRoughnessCorrection:
                        roughnessCorrectionFactor = self.roughnessCorrectionFactor(coolantReynoldsNumber, coolantPrandtlNumber, roughSmoothRatio)
                    else:
                        roughnessCorrectionFactor = 1

                    coolantNusseltNumber = self.superCriticalEquation(coolantState.cp, coolantState.viscosity, coolantPrandtlNumber, massFlowRate, self.coolingChannels.channelInstance.hydraulicDiameter, coolantState.T, coolantSideWallTemp)

                    coolantSideHeatTransferCoefficient = self.nusseltNumberToHeatTransferCoefficient(coolantNusseltNumber, self.coolingChannels.channelInstance.midWidth, coolantState.thermalConductivity, curvatureCorrection=curvatureCorrectionFactor, roughnessCorrection=roughnessCorrectionFactor)

                    # Calculate fin correction factor if included
                    if self.includeFinCorrection:
                        finCorrectionFactor = self.finCorrectionFactor(coolantSideHeatTransferCoefficient, self.coolingChannels.midRibThickness, self.coolingChannels.wallConductivity, self.coolingChannels.channelHeight)
                        coolantSideHeatTransferCoefficient = self.finCorrection(coolantSideHeatTransferCoefficient, finCorrectionFactor, self.coolingChannels.channelInstance.midWidth, self.coolingChannels.channelHeight, self.coolingChannels.midRibThickness)

                    # Calculate the coolant side wall temperature using the heat flux calculated from the gas side
                    coolantSideWallTemp = (gasSideHeatFlux / coolantSideHeatTransferCoefficient) + coolantState.T
                    surfaceCoolantState.defineState("T", coolantSideWallTemp, "P", coolantState.P)

                elif self.coolantSideHeatTransferCorrelation.lower() == "dittus-boelter" or self.coolantSideHeatTransferCorrelation.lower() == "dittus boelter":

                    # Calculate curvature correction factor if included
                    if self.includeCurvatureCorrection:
                        curvatureCorrectionFactor = self.getCurvatureCorrectionFactor(i, coolantReynoldsNumber, self.coolingChannels.channelInstance.hydraulicDiameter)
                    else:
                        curvatureCorrectionFactor = 1

                    # Calculate roughness correction factor if included
                    if self.includeRoughnessCorrection:
                        roughnessCorrectionFactor = self.roughnessCorrectionFactor(coolantReynoldsNumber, coolantPrandtlNumber, roughSmoothRatio)
                    else:
                        roughnessCorrectionFactor = 1

                    # Calculate heat transfer coefficient for the coolant side
                    coolantNusseltNumber = self.dittusBoelterEquation(coolantReynoldsNumber, coolantPrandtlNumber)
                    coolantSideHeatTransferCoefficient = self.nusseltNumberToHeatTransferCoefficient(coolantNusseltNumber, self.coolingChannels.channelInstance.midWidth, coolantState.thermalConductivity, curvatureCorrection=curvatureCorrectionFactor, roughnessCorrection=roughnessCorrectionFactor)

                    # Calculate fin correction factor if included
                    if self.includeFinCorrection:
                        finCorrectionFactor = self.finCorrectionFactor(coolantSideHeatTransferCoefficient, self.coolingChannels.midRibThickness, self.coolingChannels.wallConductivity, self.coolingChannels.channelHeight)
                        coolantSideHeatTransferCoefficient = self.applyFinCorrection(coolantSideHeatTransferCoefficient, finCorrectionFactor, self.coolingChannels.channelInstance.midWidth, self.coolingChannels.channelHeight, self.coolingChannels.midRibThickness)
                    
                    # Calculate the coolant side wall temperature using the heat flux calculated from the gas side
                    coolantSideWallTemp = (gasSideHeatFlux / coolantSideHeatTransferCoefficient) + coolantState.T
                    surfaceCoolantState.defineState("T", coolantSideWallTemp, "P", coolantState.P)
            
                elif self.coolantSideHeatTransferCorrelation.lower() == "sieder-tate" or self.coolantSideHeatTransferCorrelation.lower() == "sieder tate":

                    # Initialise the state of the coolant at the station using the value of the previous station
                    coolantSideWallTemp = coolantSideWallTempPrev
                    surfaceCoolantState.defineState("T", coolantSideWallTemp, "P", coolantState.P)

                    # Reset previous coolant side wall temperature to zero to avoid breaking the loop immediately
                    coolantSideWallTempPrev = 0

                    # Iterate to find the coolant side wall temperature as the surface viscosity is dependent on the surface temperature
                    while abs(coolantSideWallTemp - coolantSideWallTempPrev) > convergenceCriteria:

                        coolantSideWallTempPrev = coolantSideWallTemp

                        # Calculate curvature correction factor if included
                        if self.includeCurvatureCorrection:
                            curvatureCorrectionFactor = self.getCurvatureCorrectionFactor(i, coolantReynoldsNumber, self.coolingChannels.channelInstance.hydraulicDiameter)
                        else:
                            curvatureCorrectionFactor = 1

                        # Calculate roughness correction factor if included
                        if self.includeRoughnessCorrection:
                            roughnessCorrectionFactor = self.roughnessCorrectionFactor(coolantReynoldsNumber, coolantPrandtlNumber, roughSmoothRatio)
                        else:
                            roughnessCorrectionFactor = 1

                        # Calculate heat transfer coefficient for the coolant side
                        coolantNusseltNumber = self.siederTateEquation(coolantReynoldsNumber, coolantPrandtlNumber, coolantState.viscosity, surfaceCoolantState.viscosity)

                        coolantSideHeatTransferCoefficient = self.nusseltNumberToHeatTransferCoefficient(coolantNusseltNumber, self.coolingChannels.channelInstance.hydraulicDiameter, coolantState.thermalConductivity, curvatureCorrection=curvatureCorrectionFactor, roughnessCorrection=roughnessCorrectionFactor)
                        
                        # Calculate fin correction factor if included
                        if self.includeFinCorrection:
                            finCorrectionFactor = self.finCorrectionFactor(coolantSideHeatTransferCoefficient, self.coolingChannels.midRibThickness, self.coolingChannels.wallConductivity, self.coolingChannels.channelHeight)
                            coolantSideHeatTransferCoefficient = self.applyFinCorrection(coolantSideHeatTransferCoefficient, finCorrectionFactor, self.coolingChannels.channelInstance.midWidth, self.coolingChannels.channelHeight, self.coolingChannels.midRibThickness)

                        # Calculate the coolant side wall temperature using the heat flux calculated from the gas side
                        coolantSideWallTemp = (gasSideHeatFlux / coolantSideHeatTransferCoefficient) + coolantState.T
                        surfaceCoolantState.defineState("T", coolantSideWallTemp, "P", coolantState.P)

                elif self.coolantSideHeatTransferCorrelation.lower() == "ruan-meng" or self.coolantSideHeatTransferCorrelation.lower() == "ruan meng":

                    # Initialise the state of the coolant at the station using the value of the previous station
                    coolantSideWallTemp = coolantSideWallTempPrev
                    surfaceCoolantState.defineState("T", coolantSideWallTemp, "P", coolantState.P)

                    # Reset previous coolant side wall temperature to zero to avoid breaking the loop immediately
                    coolantSideWallTempPrev = 0

                    # Iterate to find the coolant side wall temperature as the surface viscosity is dependent on the surface temperature
                    while abs(coolantSideWallTemp - coolantSideWallTempPrev) > convergenceCriteria:

                        coolantSideWallTempPrev = coolantSideWallTemp

                        # Calculate curvature correction factor if included
                        if self.includeCurvatureCorrection:
                            curvatureCorrectionFactor = self.getCurvatureCorrectionFactor(i, coolantReynoldsNumber, self.coolingChannels.channelInstance.hydraulicDiameter)
                        else:
                            curvatureCorrectionFactor = 1

                        # Calculate roughness correction factor if included
                        if self.includeRoughnessCorrection:
                            roughnessCorrectionFactor = self.roughnessCorrectionFactor(coolantReynoldsNumber, coolantPrandtlNumber, roughSmoothRatio)
                        else:
                            roughnessCorrectionFactor = 1

                        # Calculate heat transfer coefficient for the coolant side
                        distance = self.thrustChamber.axialCoords[0] - self.thrustChamber.axialCoords[i]
                        coolantNusseltNumber = self.ruanMengEquation(coolantReynoldsNumber, coolantPrandtlNumber, coolantState.D, surfaceCoolantState.D, self.coolingChannels.channelInstance.hydraulicDiameter, distance)

                        coolantSideHeatTransferCoefficient = self.nusseltNumberToHeatTransferCoefficient(coolantNusseltNumber, self.coolingChannels.channelInstance.midWidth, coolantState.thermalConductivity, curvatureCorrection=curvatureCorrectionFactor, roughnessCorrection=roughnessCorrectionFactor)

                        # Calculate fin correction factor if included
                        if self.includeFinCorrection:
                            finCorrectionFactor = self.finCorrectionFactor(coolantSideHeatTransferCoefficient, self.coolingChannels.midRibThickness, self.coolingChannels.wallConductivity, self.coolingChannels.channelHeight)
                            coolantSideHeatTransferCoefficient = self.applyFinCorrection(coolantSideHeatTransferCoefficient, finCorrectionFactor, self.coolingChannels.channelInstance.midWidth, self.coolingChannels.channelHeight, self.coolingChannels.midRibThickness)

                        # Calculate the coolant side wall temperature using the heat flux calculated from the gas side
                        coolantSideWallTemp = (gasSideHeatFlux / coolantSideHeatTransferCoefficient) + coolantState.T
                        surfaceCoolantState.defineState("T", coolantSideWallTemp, "P", coolantState.P)

                # Set coolantSideWallTempPrev to the coverged value so that it can be used to initialise the next station
                # for decreased solve time
                coolantSideWallTempPrev = coolantSideWallTemp

                # Calculate the gas side wall temperature using the gas side heat flux and the coolant side wall temperature
                newGasSideWallTemp = (gasSideHeatFlux * (self.coolingChannels.wallThickness / self.coolingChannels.wallConductivity)) + coolantSideWallTemp

                # Set the gas side wall temperature guess value to the average of the previous guess value, and the calculated value above
                gasSideWallTemp = (gasSideWallTemp + newGasSideWallTemp) / 2

                print(gasSideWallTemp)

                iterations += 1

            #stationLength = distanceBetweenTwoPoints([self.thrustChamber.axialCoords[i], self.thrustChamber.radialCoords[i]],
                                                    #[self.thrustChamber.axialCoords[i-1], self.thrustChamber.radialCoords[i-1]])

            stationLength = self.thrustChamber.axialCoords[i] - self.thrustChamber.axialCoords[i-1]

            pressureLoss = self.pressureLoss(frictionFactor, stationLength, self.coolingChannels.channelInstance.hydraulicDiameter, coolantState.D, coolantVelocity)

            print("Solved station " + str(i) + " in " + str(iterations) + " iterations")
            if self.thrustChamber.axialCoords[i] == 0:
                print("\n*** Throat ***\n")
            print("Coolant Pressure: " + str(coolantState.P/1e5))
            print("Coolant Density: " + str(coolantState.D))
            print("Coolant Velocity: " + str(coolantVelocity))
            print("Adiabatic Wall Temp: " + str(adiabaticWallTemp))
            print("Gas Side Wall Temp: " + str(gasSideWallTemp))
            print("Coolant Side Wall Temp: " + str(coolantSideWallTemp))
            print("Coolant Temp: " + str(coolantState.T))
            print("Gas Side H: " + str(gasSideHeatTransferCoefficient))
            print("Coolant Side H: " + str(coolantSideHeatTransferCoefficient))
            print("Heat Flux: "+ str(gasSideHeatFlux/1e6))
                
            self.heatFluxes[i] = gasSideHeatFlux
            self.adiabaticWallTemps[i] = adiabaticWallTemp
            self.gasSideWallTemps[i] = gasSideWallTemp
            self.coolantSideWallTemps[i] = coolantSideWallTemp
            self.coolantBulkTemps[i] = coolantState.T
            self.coolantPressures[i] = coolantState.P
            self.coolantReynoldsNumbers[i] = coolantReynoldsNumber
            self.coolantNusseltNumbers[i] = coolantNusseltNumber
            self.coolantPrandtlNumbers[i] = coolantPrandtlNumber
            self.gasSideHeatTransferCoefficients[i] = gasSideHeatTransferCoefficient
            self.coolantSideHeatTransferCoefficients[i] = coolantSideHeatTransferCoefficient
            self.channelAreas[i] = self.coolingChannels.channelInstance.individualChannelArea


            temperatureChange = ((gasSideHeatFlux * 2 * pi * self.thrustChamber.radialCoords[i] * stationLength) / (self.massFlowRate * coolantState.cp))
            newCoolantBulkTemp = coolantState.T + temperatureChange
            
            newCoolantPressure = coolantState.P - pressureLoss

            coolantState.defineState("T", newCoolantBulkTemp, "P", newCoolantPressure)

            totalIterations += iterations

            i -= 1

        self.outletState.defineState("T", self.coolantBulkTemps[1], "P", self.coolantPressures[1])
        self.enthalpyChange = self.outletState.H - self.inletState.H
        self.totalHeatPower = self.massFlowRate * self.enthalpyChange

        print("---")
        print("Heat transfer calculaton complete")
        print("Total iterations: " + str(totalIterations))
        print("---")
        print("Coolant Enthalpy Change: " + str(self.enthalpyChange))
        print("Total Heat Power: " + str(self.totalHeatPower))
        print("Total Pressure Loss: " + str(self.outletState.P - self.inletState.P) + " Pa")

    def bartzEquation(self, localArea, gasSideWallTemp, localMachNumber, section):

        if section == "nozzle":

            viscosity = self.getLocalGasViscosity(localArea)
            prandtlNumber = self.getLocalGasPrandtlNumber(localArea)
            specificHeat = self.getLocalGasSpecificHeat(localArea)
            gamma = self.getLocalGamma(localArea)

        else:
            
            viscosity = self.thrustChamber.chamberViscosity
            prandtlNumber = self.thrustChamber.chamberPrandtlNumber
            specificHeat = self.thrustChamber.chamberHeatCapacity
            gamma = self.thrustChamber.chamberGamma

        firstBracket = 0.026 / (self.thrustChamber.throatDiameter ** 0.2)
        secondBracket = (((viscosity ** 0.2) * specificHeat) / (prandtlNumber ** 0.6))
        thirdBracket = ((self.thrustChamber.chamberPressureSI) / self.thrustChamber.cStar) ** 0.8
        fourthBracket = (self.thrustChamber.throatDiameter / self.thrustChamber.throatAverageRadiusOfCurvature) ** 0.1
        fifthBracket = (self.thrustChamber.throatArea / localArea) ** 0.9
        sigma = self.bartzCorrectionFactor(gasSideWallTemp, localMachNumber, gamma)

        return firstBracket * secondBracket * thirdBracket * fourthBracket * fifthBracket * sigma

    def bartzCorrectionFactor(self, gasSideWallTemp, localMachNumber, localGamma):

        denominatorFirstBracket = ((0.5 * (gasSideWallTemp / self.thrustChamber.chamberTemp) * (1 + ((localGamma - 1) / 2) * (localMachNumber ** 2)) + 0.5) ** 0.68)
        denominatorSecondBracket = (1 + ((localGamma - 1) / 2) * (localMachNumber ** 2)) ** 0.12

        return (1 / (denominatorFirstBracket * denominatorSecondBracket))

    def colebrookEquation(self, channelRoughness, hydraulicDiameter, reynoldsNumber, convergenceCriteria=0.01):

        x = -2 * log10(channelRoughness / (3.7 * hydraulicDiameter))
        xPrev = x + 10

        while abs(abs(x) - abs(xPrev)) > convergenceCriteria:

            xPrev = x
            
            x = -2 * log10((channelRoughness / (3.7 * hydraulicDiameter)) + ((2.51 * xPrev) / reynoldsNumber))

        frictionFactor = (1 / x) ** 2

        return frictionFactor
        
    def pressureLoss(self, frictionFactor, length, hydraulicDiameter, density, velocity):

        return frictionFactor * (length / hydraulicDiameter) * 0.5 * density * (velocity ** 2)
    
    def reynoldsNumber(self, density, velocity, characteristicLength, visocity):

        return (density * velocity * characteristicLength) / visocity

    def prandtlNumber(self, specificHeat, viscosity, thermalConductivity):

        return (specificHeat * viscosity) / thermalConductivity

    def massFlowRateToVelocity(self, massFlowRate, density, area):

        return massFlowRate / (density * area)

    def dittusBoelterEquation(self, reynoldsNumber, prandtlNumber, C1=0.026):

        return C1 * (reynoldsNumber ** 0.8) * (prandtlNumber ** 0.4)

    def siederTateEquation(self, reynoldsNumber, prandtlNumber, bulkViscosity, surfaceViscosity, C1=0.026):

        return C1 * (reynoldsNumber ** (4/5)) * (prandtlNumber ** (1/3)) * ((bulkViscosity / surfaceViscosity) ** 0.14)

    def superCriticalEquation(self, specificHeat, viscosity, prandtlNumber, massFlowRate, hydraulicDiameter, coolantBulkTemp, coolantSideWallTemp):

        return ((0.026 * specificHeat * (viscosity ** 0.2)) / (prandtlNumber ** (2/3))) * (((massFlowRate * G) ** 0.8) / (hydraulicDiameter ** 0.2)) * ((coolantBulkTemp / coolantSideWallTemp) ** 0.55)

    def taylorEquation(self, reynoldsNumber, prandtlNumber, coolantBulkTemp, coolantSurfaceTemp, hydraulicDiameter, distance):

        exponent = 0.57 - 1.59 * (hydraulicDiameter / distance)

        return 0.023 * (reynoldsNumber ** 0.8) * (prandtlNumber ** 0.4) * ((coolantBulkTemp / coolantSurfaceTemp) ** exponent)

    def ruanMengEquation(self, reynoldsNumber, prandtlNumber, coolantBulkDensity, coolantSurfaceDensity, hydraulicDiameter, distance):

        return 0.0069 * (reynoldsNumber ** 0.8) * (prandtlNumber ** 0.66) * ((coolantBulkDensity / coolantSurfaceDensity) ** 0.43) * (1 + 2.4 * (hydraulicDiameter / distance))

    def nusseltNumberToHeatTransferCoefficient(self, nusseltNumber, charactersiticLength, thermalConductivity, curvatureCorrection=1, roughnessCorrection=1):

        return (nusseltNumber * thermalConductivity * curvatureCorrection * roughnessCorrection) / charactersiticLength

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

        result = CEA.get_Temperatures(Pc=self.thrustChamber.injectionPressure, MR=self.thrustChamber.mixtureRatio, eps=expansionRatioAtArea, frozen=c, frozenAtThroat=tc)

        return result[2]
        
    def getLocalGasTempChamber(self, area):

        mach = self.getLocalMachNumberChamber(area)
        gamma = self.getLocalGamma(area)

        return self.thrustChamber.chamberTemp / (1 + ((gamma - 1) / 2) * (mach ** 2))

    def getLocalGasViscosity(self, area):

        expansionRatioAtArea = area / self.thrustChamber.throatArea

        result = self.thrustChamber.getExitTransportPropertiesAtExpansionRatio(expansionRatioAtArea)

        return result[1] * 0.001

    def getLocalGasPrandtlNumber(self, area):

        expansionRatioAtArea = area / self.thrustChamber.throatArea

        result = self.thrustChamber.getExitTransportPropertiesAtExpansionRatio(expansionRatioAtArea)

        return result[3]
    
    def getLocalGasSpecificHeat(self, area):

        expansionRatioAtArea = area / self.thrustChamber.throatArea

        result = self.thrustChamber.getExitTransportPropertiesAtExpansionRatio(expansionRatioAtArea)

        return result[0]

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

        return gasTemp

    def getLocalAdiabaticWallTempChamber(self, area):

        mach = self.getLocalMachNumberChamber(area)
        r = self.thrustChamber.chamberPrandtlNumber ** 0.33
        gamma = self.thrustChamber.chamberGamma
        gasTemp = self.getLocalGasTempChamber(area)

        numerator = 1 + (r * ((gamma-1) / 2) * (mach ** 2))
        denominator = 1 + ((gamma-1) / 2) * (mach ** 2)

        return gasTemp

    def finCorrectionFactor(self, heatTransferCoeffcient, ribThickness, thermalConductivity, channelHeight):
        
        temp = sqrt((2 * heatTransferCoeffcient * ribThickness) / thermalConductivity) * (channelHeight / ribThickness)

        return tanh(temp) / temp

    def applyFinCorrection(self, heatTransferCoefficient, finCorrectionFactor, channelWidth, channelHeight, ribThickness):

        numerator = channelWidth + 2 * finCorrectionFactor * channelHeight
        denominator = channelWidth + ribThickness

        return heatTransferCoefficient * (numerator / denominator)

    def curvatureCorrection(self, p1, p2, p3, reynoldsNumber, hydraulicDiameter):

        radiusOfCurvature, centrePoint = radiusOfCurvature3Points2D(p1, p2, p3, returnCenterPoint=True)

        if radiusOfCurvature == None:

            return 1

        if centrePoint[1] > p2[1]:
            exponent = 0.05
        else:
            exponent = -0.05

        return (reynoldsNumber * ((hydraulicDiameter / (2 * radiusOfCurvature)) ** 2)) ** exponent

    def smoothFrictionFactor(self, reynoldsNumber):

        if reynoldsNumber >= 10 ** 4:

            return 0.0014 + (0.125 / (reynoldsNumber ** 0.32))

        elif reynoldsNumber < 10 ** 4:

            return 0.0791 / (reynoldsNumber ** 0.25)

    def roughnessCorrectionFactor(self, reynoldsNumber, prandtlNumber, roughSmoothRatio):

        numerator = 1 + 1.5 * (prandtlNumber ** (-1/6)) * (reynoldsNumber ** (-1/8)) * (prandtlNumber - 1)
        denominator = 1 + 1.5 * (prandtlNumber ** (-1/6)) * (reynoldsNumber ** (-1/8)) * ((prandtlNumber * roughSmoothRatio) - 1)

        return (numerator / denominator) * roughSmoothRatio

    def getCurvatureCorrectionFactor(self, station, reynoldsNumber, hydraulicDiameter, included=True):

        prevCoord = [self.thrustChamber.axialCoords[station-1], self.thrustChamber.radialCoords[station-1]]
        currentCoord = [self.thrustChamber.axialCoords[station], self.thrustChamber.radialCoords[station]]
        nextCoord = [self.thrustChamber.axialCoords[station+1], self.thrustChamber.radialCoords[station+1]]

        return self.curvatureCorrection(prevCoord, currentCoord, nextCoord, reynoldsNumber, hydraulicDiameter)
