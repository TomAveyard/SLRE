import sys
from PropTools.SubSystems.Engine.Propellant.propellant import Propellant
import numpy as np
from math import pi, sqrt, sin, radians
from PropTools.SubSystems.Engine.ThrustChamber.thrustChamber import ThrustChamber
from PropTools.Utils.mathsUtils import radiusToArea, distanceBetweenTwoPoints
from PropTools.SubSystems.Engine.Cycle.component import Component
from PropTools.Thermo.dimensionlessNumbers import reynoldsNumber, prandtlNumber
import PropTools.Thermo.heatTransfer as ht
import PropTools.Thermo.fluidDynamics as fd
import matplotlib.pyplot as plt
from PropTools.SubSystems.Engine.ThrustChamber.coolingChannels import CoolingChannels
from PropTools.Utils.helixGeometry import Helix


class SolverParameters:

    def __init__(self, 
    bartzEquationCoefficient: float = 0.023,
    coolantSideHeatTransferCorrelation: str = "dittus-boelter",
    coolantSideHeatTransferCorrelationCoefficient: str = "default",
    includeFinCorrection: bool = True,
    includeCurvatureCorrection: bool = True,
    includeRoughnessCorrection: bool = True,
    convergenceCriteria: float = 0.01):

        self.bartzEquationCoefficient = bartzEquationCoefficient
        self.coolantSideHeatTransferCorrelation = coolantSideHeatTransferCorrelation
        self.coolantSideHeatTransferCorrelationCoefficient = coolantSideHeatTransferCorrelationCoefficient
        self.includeFinCorrection = includeFinCorrection
        self.includeCurvatureCorrection = includeCurvatureCorrection
        self.includeRoughnessCorrection = includeRoughnessCorrection
        self.convergenceCriteria = convergenceCriteria


# Class to store information about the regenerative cooling system
# Performs an iterative calculation for every point defined by the input thrust chamber object
# with the channel dimensions defined by the cooling channel object
class RegenerativeCooling(Component):

    def __init__(self, thrustChamber: ThrustChamber = None, coolingChannels: CoolingChannels = None, solverParameters: SolverParameters = None):

        super().__init__()
        self.type = "regenerative cooling"

        self.thrustChamber = thrustChamber
        self.coolingChannels = coolingChannels
        self.solverParameters = solverParameters

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

    def calculate(self, inletState: Propellant = None, massFlowRate: float = None) -> None:

        print("---")
        print("Starting Heat Transfer Calculation")
        print("Number Of Stations: " + str(self.numberOfStations))

        self.inletState = inletState
        self.massFlowRate = massFlowRate

        # Object for the overall outlet state
        self.outletState = Propellant(self.inletState.name)

        # Get CEA object that we can query
        cea = self.thrustChamber.getCEAObject()
        c, ct = self.thrustChamber.convertConditionToCEAFlag(self.thrustChamber.condition, self.thrustChamber.throatCondition)
        injectionPressure = self.thrustChamber.injectionPressure
        mixtureRatio = self.thrustChamber.mixtureRatio
        exitExpansionRatio = self.thrustChamber.expansionRatio

        # Get chamber properties
        chamberTemp = self.thrustChamber.chamberTemp
        chamberPressure = self.thrustChamber.chamberPressure * 1e5
        chamberGamma = cea.get_Chamber_MolWt_gamma(Pc=injectionPressure, MR=mixtureRatio, eps=exitExpansionRatio)[1]
        chamberSpecificHeat = cea.get_HeatCapacities(Pc=injectionPressure, MR=mixtureRatio, eps=exitExpansionRatio, frozen=c, frozenAtThroat=ct)[0]
        chamberTransport = cea.get_Chamber_Transport(Pc=injectionPressure, MR=mixtureRatio, eps=exitExpansionRatio, frozen=c)
        chamberViscosity = chamberTransport[1] * 1e-3
        chamberThermalConductivity = chamberTransport[2] * 100
        chamberPrandtlNumber = chamberTransport[3]
        cStar = cea.get_Cstar(Pc=injectionPressure, MR=mixtureRatio)

        totalIterations = 0

        # -1 to account for 0 index
        station = self.numberOfStations - 1

        # State trackers for the coolant at the inlet, outlet, and surface of a station
        stationInletState = Propellant(self.inletState.name)
        stationOutletState = Propellant(self.inletState.name)
        stationSurfaceState = Propellant(self.inletState.name)
        stationInletState.defineState("T", self.inletState.T, "P", self.inletState.P)
        stationOutletState.defineState("T", self.inletState.T, "P", self.inletState.P)

        # Counter for total distance coolant has travelled, required for some heat transfer correlations
        distance = 0
        
        # Iterate through each station from end of nozzle to start of chamber
        while station > 0:

            loopIterations = 0

            # Get coords of station
            stationAxialCoord = self.thrustChamber.axialCoords[station]
            stationRadialCoord = self.thrustChamber.radialCoords[station]

            # Flag whether we are in the nozzle or chamber
            if stationAxialCoord >= 0:
                nozzle = True
            else:
                nozzle = False

            # Get thrust chamber dimensions at the station
            stationRadius = stationRadialCoord
            stationDiameter = stationRadius * 2
            stationArea = radiusToArea(stationRadius)

            # Get throat dimensions
            throatRadius = self.thrustChamber.throatRadius
            throatDiameter = throatRadius * 2
            throatArea = radiusToArea(throatRadius)

            stationExpansionRatio = stationArea / throatArea

            # Update channel dimensions
            self.coolingChannels.channelInstance.getChannelDimensions(stationArea)

            # Get the gas properties at the station
            if nozzle:
                gasGamma = cea.get_exit_MolWt_gamma(Pc=injectionPressure, MR=mixtureRatio, eps=stationExpansionRatio)[1]
                gasMachNumber = cea.get_MachNumber(Pc=injectionPressure, MR=mixtureRatio, eps=stationExpansionRatio, frozen=c, frozenAtThroat=ct)
                gasTemp = cea.get_Temperatures(Pc=injectionPressure, MR=mixtureRatio, eps=stationExpansionRatio, frozen=c, frozenAtThroat=ct)[-1]
                gasSpecificHeat = cea.get_HeatCapacities(Pc=injectionPressure, MR=mixtureRatio, eps=stationExpansionRatio, frozen=c, frozenAtThroat=ct)[-1]
                gasTransport = cea.get_Exit_Transport(Pc=injectionPressure, MR=mixtureRatio, eps=stationExpansionRatio, frozen=c)
                gasViscosity = gasTransport[1] * 1e-3
                gasThermalConductivity = gasTransport[2] * 100
                gasPrandtlNumber = gasTransport[3]
                gasAdiabaticWallTemp = fd.getAdiabaticWallTemp(gasTemp, gasMachNumber, gasGamma, gasPrandtlNumber)
            else:
                gasGamma = cea.get_Chamber_MolWt_gamma(Pc=injectionPressure, MR=mixtureRatio, eps=stationExpansionRatio)[1]
                gasMachNumber = fd.getSubsonicMachNumber(stationArea, throatArea, gasGamma, accuracy=self.solverParameters.convergenceCriteria)
                gasTemp = fd.getGasTemp(stationArea, gasMachNumber, gasGamma, chamberTemp)
                gasSpecificHeat = cea.get_HeatCapacities(Pc=injectionPressure, MR=mixtureRatio, eps=stationExpansionRatio, frozen=c, frozenAtThroat=ct)[0]
                gasTransport = cea.get_Chamber_Transport(Pc=injectionPressure, MR=mixtureRatio, eps=stationExpansionRatio, frozen=c)
                gasViscosity = gasTransport[1] * 1e-3
                gasThermalConductivity = gasTransport[2] * 100
                gasPrandtlNumber = gasTransport[3]
                gasAdiabaticWallTemp = fd.getAdiabaticWallTemp(gasTemp, gasMachNumber, gasGamma, gasPrandtlNumber)

            # Set the coolant state to the state of the outlet of the previous station
            stationInletState.defineState("T", stationOutletState.T, "P", stationOutletState.P)

            # Get coolant properties at the station inlet
            coolantVelocity = fd.massFlowRateToVelocity(self.massFlowRate, stationInletState.D, self.coolingChannels.channelInstance.totalChannelArea)
            coolantReynoldsNumber = reynoldsNumber(stationInletState.D, coolantVelocity, self.coolingChannels.channelInstance.hydraulicDiameter, stationInletState.viscosity)
            coolantPrandtlNumber = prandtlNumber(stationInletState.cp, stationInletState.viscosity, stationInletState.thermalConductivity)

            # Get friction factor
            if self.coolingChannels.wallRoughnessHeight > 0:
                frictionFactor = fd.colebrookEquation(self.coolingChannels.wallRoughnessHeight, self.coolingChannels.channelInstance.hydraulicDiameter, coolantReynoldsNumber, convergenceCriteria=self.solverParameters.convergenceCriteria)
            elif self.coolingChannels.wallRoughnessHeight == 0:
                frictionFactor = fd.smoothFrictionFactor(coolantReynoldsNumber)
            else:
                exit("Please input a positive channel roughness")

            # Initialisation for iteration loop. After the first station, the previous station is used to initialise
            if totalIterations != 0:
                gasSideWallTemp = self.gasSideWallTemps[station+1]
                coolantSideWallTemp = self.coolantSideWallTemps[station+1]
            else:
                gasSideWallTemp = gasAdiabaticWallTemp * 0.8
                coolantSideWallTemp = stationInletState.T

            # Change the new temps so that the loop doesn't break
            newGasSideWallTemp = gasSideWallTemp + 10
            newCoolantSideWallTemp = coolantSideWallTemp + 10

            while abs(abs(gasSideWallTemp) - abs(newGasSideWallTemp)) > self.solverParameters.convergenceCriteria and abs(abs(coolantSideWallTemp) - abs(newCoolantSideWallTemp)) > self.solverParameters.convergenceCriteria:

                gasSideWallTemp = newGasSideWallTemp
                coolantSideWallTemp = newCoolantSideWallTemp

                gasSideHeatTransferCoefficient = ht.bartzEquation(throatDiameter, chamberViscosity, chamberSpecificHeat, chamberPrandtlNumber, chamberPressure, chamberTemp, cStar, stationArea, gasSideWallTemp, gasMachNumber, gasGamma, C1=self.solverParameters.bartzEquationCoefficient)
                stationSurfaceState.defineState("T", coolantSideWallTemp, "P", stationOutletState.P)

                # Get nusselt number using chosen correlation
                if self.solverParameters.coolantSideHeatTransferCorrelation == "dittus-boelter" or self.solverParameters.coolantSideHeatTransferCorrelation == "dittus boelter":

                    if self.solverParameters.coolantSideHeatTransferCorrelationCoefficient == "default":
                        coolantNusseltNumber = ht.dittusBoelterEquation(coolantReynoldsNumber, coolantPrandtlNumber)
                    else:
                        coolantNusseltNumber = ht.dittusBoelterEquation(coolantReynoldsNumber, coolantPrandtlNumber, C1=self.solverParameters.coolantSideHeatTransferCorrelationCoefficient)

                elif self.solverParameters.coolantSideHeatTransferCorrelation == "sieder-tate" or self.solverParameters.coolantSideHeatTransferCorrelation == "sieder tate":
                    
                    if self.solverParameters.coolantSideHeatTransferCorrelationCoefficient == "default":
                        coolantNusseltNumber = ht.siederTateEquation(coolantReynoldsNumber, coolantPrandtlNumber, stationInletState.viscosity, stationSurfaceState.viscosity)
                    else:
                        coolantNusseltNumber = ht.siederTateEquation(coolantReynoldsNumber, coolantPrandtlNumber, stationInletState.viscosity, stationSurfaceState.viscosity, C1=self.solverParameters.coolantSideHeatTransferCorrelationCoefficient)

                elif self.solverParameters.coolantSideHeatTransferCorrelation == "taylor":

                    if self.solverParameters.coolantSideHeatTransferCorrelationCoefficient == "default":
                        coolantNusseltNumber = ht.taylorEquation(coolantReynoldsNumber, coolantPrandtlNumber, stationInletState.T, stationSurfaceState.T, self.coolingChannels.channelInstance.hydraulicDiameter, distance)
                    else:
                        coolantNusseltNumber = ht.siederTateEquation(coolantReynoldsNumber, coolantPrandtlNumber, stationInletState.viscosity, stationSurfaceState.viscosity, C1=self.solverParameters.coolantSideHeatTransferCorrelationCoefficient)

                elif self.solverParameters.coolantSideHeatTransferCorrelation == "ruan-meng" or self.solverParameters.coolantSideHeatTransferCorrelation == "ruan meng":
                    
                    if self.solverParameters.coolantSideHeatTransferCorrelationCoefficient == "default":
                        coolantNusseltNumber = ht.ruanMengEquation(coolantReynoldsNumber, coolantPrandtlNumber, stationInletState.D, stationSurfaceState.D, self.coolingChannels.channelInstance.hydraulicDiameter, distance)
                    else:
                        coolantNusseltNumber = ht.siederTateEquation(coolantReynoldsNumber, coolantPrandtlNumber, stationInletState.viscosity, stationSurfaceState.viscosity, C1=self.solverParameters.coolantSideHeatTransferCorrelationCoefficient)
                    
                else:
                    exit("Please choose a valid heat transfer correlation for the coolant")

                # Apply curvature correction if included
                if self.solverParameters.includeCurvatureCorrection:
                    # This will fail on the first station as there is no point before it, so 1 is returned in that case
                    try:
                        curvatureCorrectionFactor = ht.curvatureCorrectionFactor(self.thrustChamber.axialCoords[station+1], stationAxialCoord, self.thrustChamber.axialCoords[station-1], coolantReynoldsNumber, self.coolingChannels.channelInstance.hydraulicDiameter)
                    except:
                        curvatureCorrectionFactor = 1

                    coolantNusseltNumber = coolantNusseltNumber * curvatureCorrectionFactor

                # Apply roughness correction if included
                if self.solverParameters.includeRoughnessCorrection:

                    roughnessCorrectionFactor = ht.roughnessCorrectionFactor(frictionFactor, coolantReynoldsNumber, coolantPrandtlNumber)
                    coolantNusseltNumber = coolantNusseltNumber * roughnessCorrectionFactor

                # Calculate the coolant side heat transfer coefficient
                coolantSideHeatTransferCoefficient = ht.nusseltNumberToHeatTransferCoefficient(coolantNusseltNumber, self.coolingChannels.channelInstance.hydraulicDiameter, stationInletState.thermalConductivity)

                # Apply fin correction if included
                if self.solverParameters.includeFinCorrection:

                    finEffectiveness = ht.finEffectiveness(coolantSideHeatTransferCoefficient, self.coolingChannels.ribThickness, self.coolingChannels.wallConductivity, self.coolingChannels.channelHeight)
                    coolantSideHeatTransferCoefficient = ht.applyFinEffectiveness(coolantSideHeatTransferCoefficient, finEffectiveness, self.coolingChannels.channelInstance.midWidth, self.coolingChannels.channelHeight, self.coolingChannels.ribThickness)

                heatFlux = (gasAdiabaticWallTemp - stationInletState.T) / ((1 / gasSideHeatTransferCoefficient) + (self.coolingChannels.wallThickness / self.coolingChannels.wallConductivity) + (1 / coolantSideHeatTransferCoefficient))

                newGasSideWallTemp = gasAdiabaticWallTemp - (heatFlux / gasSideHeatTransferCoefficient)
                newCoolantSideWallTemp = stationInletState.T + (heatFlux / coolantSideHeatTransferCoefficient)

                loopIterations += 1

            gasSideWallTemp = newGasSideWallTemp
            coolantSideWallTemp = newCoolantSideWallTemp

            # Calculate length of the station
            stationLength = distanceBetweenTwoPoints([stationAxialCoord, stationRadialCoord], [self.thrustChamber.axialCoords[station-1], self.thrustChamber.radialCoords[station-1]]) / sin(radians(self.coolingChannels.helixAngle))
            distance += stationLength

            # Calculate temperature change over the station
            stationHeatTransferArea = 2 * pi * stationRadius * stationLength
            stationTempChange = (heatFlux * stationHeatTransferArea) / (massFlowRate * stationInletState.cp)

            # Calculate pressure loss over the station
            stationPressureLoss = fd.pressureLoss(frictionFactor, stationLength, self.coolingChannels.channelInstance.hydraulicDiameter, stationInletState.D, coolantVelocity)

            # Calculate the state of the station outlet
            stationOutletState.defineState("T", stationInletState.T + stationTempChange ,"P", stationInletState.P - stationPressureLoss)

            # Add station variables to lists
            self.heatFluxes[station] = heatFlux
            self.adiabaticWallTemps[station] = gasAdiabaticWallTemp
            self.gasSideWallTemps[station] = gasSideWallTemp
            self.coolantSideWallTemps[station] = coolantSideWallTemp
            self.coolantBulkTemps[station] = stationOutletState.T
            self.coolantPressures[station] = stationOutletState.P
            self.coolantReynoldsNumbers[station] = coolantReynoldsNumber
            self.coolantNusseltNumbers[station] = coolantNusseltNumber
            self.coolantPrandtlNumbers[station] = coolantPrandtlNumber
            self.gasSideHeatTransferCoefficients[station] = gasSideHeatTransferCoefficient
            self.coolantSideHeatTransferCoefficients[station] = coolantSideHeatTransferCoefficient
            self.channelAreas[station] = self.coolingChannels.channelInstance.individualChannelArea

            print("---")
            print("Solved station " + str(station) + " in " + str(loopIterations) + " iterations")
            print("---")
            print("Coolant Pressure: " + str(round(stationOutletState.P/1e5, 2)))
            print("Coolant Density: " + str(round(stationOutletState.D, 2)))
            print("Coolant Velocity: " + str(round(coolantVelocity, 2)))
            print("Adiabatic Wall Temp: " + str(round(gasAdiabaticWallTemp, 2)))
            print("Gas Side Wall Temp: " + str(round(gasSideWallTemp, 2)))
            print("Coolant Side Wall Temp: " + str(round(coolantSideWallTemp, 2)))
            print("Coolant Temp: " + str(round(stationOutletState.T, 2)))
            print("Gas Side H: " + str(round(gasSideHeatTransferCoefficient, 2)))
            print("Coolant Side H: " + str(round(coolantSideHeatTransferCoefficient, 2)))
            print("Heat Flux: "+ str(round(heatFlux/1e6, 2)))
            if self.solverParameters.includeCurvatureCorrection:
                print("Curvature Correction: " + str(round(curvatureCorrectionFactor, 2)))
            if self.solverParameters.includeRoughnessCorrection:
                print("Roughness Correction: " + str(round(roughnessCorrectionFactor, 2)))
            if self.solverParameters.includeFinCorrection:
                print("Fin Effectiveness: " + str(round(finEffectiveness, 2)))
            
            station -= 1
            totalIterations += loopIterations

        self.outletState.defineState("T", stationOutletState.T, "P", stationOutletState.P)
        self.enthalpyChange = self.outletState.H - self.inletState.H
        self.totalHeatPower = self.massFlowRate * self.enthalpyChange

        print("---")
        print("Heat transfer calculaton complete")
        print("Total iterations: " + str(totalIterations))
        print("---")
        print("Coolant Outlet Temperature: " + str(self.outletState.T))
        print("Coolant Enthalpy Change: " + str(self.enthalpyChange))
        print("Total Heat Power: " + str(self.totalHeatPower))
        print("Total Pressure Loss: " + str((self.outletState.P - self.inletState.P)/1e5) + " Bar")

    def plotChannels(self, plotAll: bool = True, includeThrustChamberContour: bool = False, lineWidth: float = 1):

        geometry = Helix(self.thrustChamber.axialCoords, self.thrustChamber.radialCoords, self.coolingChannels.helixAngle)
        if  plotAll:
            numberOfChannels = self.coolingChannels.numberOfChannels
        else:
            numberOfChannels = 1
        
        geometry.plotHelix(include2DContour=includeThrustChamberContour, numberOfChannels=numberOfChannels, lw=lineWidth)

    def plotHeatFlux(self, show: bool = True, save: bool = False, indexRange: list = [1, -1]) -> None:

        fig, ax = plt.subplots()
        ax.plot(self.thrustChamber.axialCoords[indexRange[0]:indexRange[1]], self.heatFluxes[indexRange[0]:indexRange[1]]/1e6)
        ax.set_xlabel("Axial Distance [m]")
        ax.set_ylabel("Heat Flux [MW m^-2]")
        ax.set_title("Thrust Chamber Heat Flux")

        if show:
            plt.show()

    def plotAdiabaticWallTemp(self, show: bool = True, save: bool = False, indexRange: list = [1, -1]) -> None:

        fig, ax = plt.subplots()
        ax.plot(self.thrustChamber.axialCoords[indexRange[0]:indexRange[1]], self.adiabaticWallTemps[indexRange[0]:indexRange[1]])
        ax.set_xlabel("Axial Distance [m]")
        ax.set_ylabel("Adiabatic Wall Temperature [K]")
        ax.set_title("Thrust Chamber Adiabatic Wall Temperature")
        
        if show:
            plt.show()

    def plotGasSideWallTemp(self, show: bool = True, save: bool = False, indexRange: list = [1, -1]) -> None:

        fig, ax = plt.subplots()
        ax.plot(self.thrustChamber.axialCoords[indexRange[0]:indexRange[1]], self.gasSideWallTemps[indexRange[0]:indexRange[1]])
        ax.set_xlabel("Axial Distance [m]")
        ax.set_ylabel("Gas Side Wall Temperature [K]")
        ax.set_title("Thrust Chamber Gas Side Wall Temperature")
        
        if show:
            plt.show()

    def plotCoolantSideWallTemp(self, show: bool = True, save: bool = False, indexRange: list = [1, -1]) -> None:

        fig, ax = plt.subplots()
        ax.plot(self.thrustChamber.axialCoords[indexRange[0]:indexRange[1]], self.coolantSideWallTemps[indexRange[0]:indexRange[1]])
        ax.set_xlabel("Axial Distance [m]")
        ax.set_ylabel("Coolant Side Wall Temperature [K]")
        ax.set_title("Thrust Chamber Coolant Side Wall Temperature")
        
        if show:
            plt.show()

    def plotCoolantBulkTemp(self, show: bool = True, save: bool = False, indexRange: list = [1, -1]) -> None:

        fig, ax = plt.subplots()
        ax.plot(self.thrustChamber.axialCoords[indexRange[0]:indexRange[1]], self.coolantBulkTemps[indexRange[0]:indexRange[1]])
        ax.set_xlabel("Axial Distance [m]")
        ax.set_ylabel("Coolant Bulk Temperature [K]")
        ax.set_title("Thrust Chamber Coolant Bulk Temperature")
        
        if show:
            plt.show()

    def plotCoolantPressure(self, show: bool = True, save: bool = False, indexRange: list = [1, -1]) -> None:

        fig, ax = plt.subplots()
        ax.plot(self.thrustChamber.axialCoords[indexRange[0]:indexRange[1]], self.coolantPressures[indexRange[0]:indexRange[1]]/1e5)
        ax.set_xlabel("Axial Distance [m]")
        ax.set_ylabel("Coolant Pressure [bar]")
        ax.set_title("Thrust Chamber Coolant Pressure")
        
        if show:
            plt.show()

    def plotReynoldsNumber(self, show: bool = True, save: bool = False, indexRange: list = [1, -1]) -> None:

        fig, ax = plt.subplots()
        ax.plot(self.thrustChamber.axialCoords[indexRange[0]:indexRange[1]], self.coolantReynoldsNumbers[indexRange[0]:indexRange[1]])
        ax.set_xlabel("Axial Distance [m]")
        ax.set_ylabel("Coolant Reynolds Number")
        ax.set_title("Thrust Chamber Coolant Reynolds Number")
        
        if show:
            plt.show()

    def plotPrandtlNumber(self, show: bool = True, save: bool = False, indexRange: list = [1, -1]) -> None:

        fig, ax = plt.subplots()
        ax.plot(self.thrustChamber.axialCoords[indexRange[0]:indexRange[1]], self.coolantPrandtlNumbers[indexRange[0]:indexRange[1]])
        ax.set_xlabel("Axial Distance [m]")
        ax.set_ylabel("Coolant Prandtl Number")
        ax.set_title("Thrust Chamber Coolant Prandtl Number")
        
        if show:
            plt.show()

    def plotGasSideHeatTransferCoefficient(self, show: bool = True, save: bool = False, indexRange: list = [1, -1]) -> None:

        fig, ax = plt.subplots()
        ax.plot(self.thrustChamber.axialCoords[indexRange[0]:indexRange[1]], self.gasSideHeatTransferCoefficients[indexRange[0]:indexRange[1]])
        ax.set_xlabel("Axial Distance [m]")
        ax.set_ylabel("Gas Side Heat Transfer Coefficient [W m^-1 K^-1]")
        ax.set_title("Thrust Chamber Gas Side Heat Transfer Coefficient")
        
        if show:
            plt.show()

    def plotCoolantSideHeatTransferCoefficient(self, show: bool = True, save: bool = False, indexRange: list = [1, -1]) -> None:

        fig, ax = plt.subplots()
        ax.plot(self.thrustChamber.axialCoords[indexRange[0]:indexRange[1]], self.coolantSideHeatTransferCoefficients[indexRange[0]:indexRange[1]])
        ax.set_xlabel("Axial Distance [m]")
        ax.set_ylabel("Coolant Side Heat Transfer Coefficient [W m^-1 K^-1]")
        ax.set_title("Thrust Chamber Coolant Side Heat Transfer Coefficient")
        
        if show:
            plt.show()
    
    def plotChannelAreas(self, show: bool = True, save: bool = False, indexRange: list = [1, -1]) -> None:

        fig, ax = plt.subplots()
        ax.plot(self.thrustChamber.axialCoords[indexRange[0]:indexRange[1]], self.channelAreas[indexRange[0]:indexRange[1]])
        ax.set_xlabel("Axial Distance [m]")
        ax.set_ylabel("Cooling Channel Area [m^-2]")
        ax.set_title("Thrust Chamber Cooling Channel Area")
        
        if show:
            plt.show()
