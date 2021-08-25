from math import sqrt
import sys
from matplotlib import pyplot as plt
from rocketcea.cea_obj_w_units import CEA_Obj
import numpy as np

from PropTools.Utils.constants import G, R
from PropTools.Utils.mathsUtils import areaToRadius, radiusToArea
from PropTools.SubSystems.Engine.Propellant.propellant import Propellant
from PropTools.SubSystems.Engine.ThrustChamber.combustionChamber import CombustionChamber
from PropTools.SubSystems.Engine.ThrustChamber.nozzle import ConicalNozzle, RaoBellNozzle

class ThrustChamber:

    def __init__(self, 
                fuelName, 
                oxName,
                thrust, 
                chamberPressure, 
                mixtureRatioOverride=False, 
                ambientPressure=1.01325, 
                fac=False,
                CR=None,
                facPlenumPressureSpecified=True,
                mixtureRatioSearchResolution=0.1, 
                mixtureRatioSearchStart=1, 
                condition='equilibrium', 
                throatCondition='equilibrium'):

        self.fuel = Propellant(fuelName)
        self.ox = Propellant(oxName)

        self.thrust = thrust

        self.chamberPressure = chamberPressure
        self.ambientPressure = ambientPressure

        self.fac = fac
        self.CR = CR

        # If a finite area combustor contraction ratio is specified and the chamber pressure given is the plenum/combustion pressure
        # then this estimates an injector face pressure for use in the calculations.
        # Otherwise the injection pressure is the chamber pressure specified.
        if self.fac == True and facPlenumPressureSpecified == True:
            self.facPlenumPressureSpecified = facPlenumPressureSpecified
            self.injectionPressureOverCombustionPressure = 1.0 + 0.54 / self.CR**2.2
            self.injectionPressure = self.chamberPressure * self.injectionPressureOverCombustionPressure
        else:
            self.injectionPressure = self.chamberPressure

        self.condition = condition
        self.throatCondition = throatCondition

        self.CEA = self.getCEAObject()

        if mixtureRatioOverride != False:

            self.mixtureRatio = mixtureRatioOverride
            self.getCEAResults()

        else:

            self.mixtureRatio = self.getMaxSpecificImpulseMixtureRatio(searchResolution=mixtureRatioSearchResolution, startSearch=mixtureRatioSearchStart)
            self.getCEAResults()

        self.getMassFlowRate()
        self.getExitVelocity()
        self.getExitSizes()
        self.getThroatSizes()

        self.combustionChamber = None
        self.nozzle = None
        self.axialCoords = None
        self.radialCoords = None

    # Returns a cea object using the unit system seen below
    def getCEAObject(self):

        # If facCR (finite area combustor contraction ratio) is not specified, then an infinite area combustor is assumed
        if self.fac == False:

            CEAObject = CEA_Obj(oxName=self.ox.ceaName, 
                fuelName=self.fuel.ceaName, 
                isp_units = 'sec',
                cstar_units = 'm/s',
                pressure_units = 'Bar',
                temperature_units = 'K',
                sonic_velocity_units = 'm/s',
                enthalpy_units = 'J/kg',
                density_units = 'kg/m^3',
                specific_heat_units = 'J/kg-K',
                viscosity_units = 'centipoise',
                thermal_cond_units = 'W/cm-degC'
                )

        elif self.fac == True:

            CEAObject = CEA_Obj(oxName=self.ox.ceaName, 
                fuelName=self.fuel.ceaName,
                fac_CR=self.CR, 
                isp_units = 'sec',
                cstar_units = 'm/s',
                pressure_units = 'Bar',
                temperature_units = 'K',
                sonic_velocity_units = 'm/s',
                enthalpy_units = 'J/kg',
                density_units = 'kg/m^3',
                specific_heat_units = 'J/kg-K',
                viscosity_units = 'centipoise',
                thermal_cond_units = 'W/cm-degC'
                )

        return CEAObject

    # Gets the most relevant CEA results and assigns them to easy to access variables in the thrustChamber object
    def getCEAResults(self):

        c, tc = self.convertConditionToCEAFlag(self.condition, self.throatCondition)
        
        self.expansionRatio = self.CEA.get_eps_at_PcOvPe(Pc=self.injectionPressure, MR=self.mixtureRatio, PcOvPe=self.chamberPressure/self.ambientPressure, frozen=c, frozenAtThroat=tc)
        self.specificImpulse = self.CEA.estimate_Ambient_Isp(Pc=self.injectionPressure, MR=self.mixtureRatio, eps=self.expansionRatio, Pamb=self.ambientPressure, frozen=c, frozenAtThroat=tc)[0]
        self.specificImpulseVacuum = self.CEA.get_Isp(Pc=self.injectionPressure, MR=self.mixtureRatio, eps=self.expansionRatio, frozen=c, frozenAtThroat=tc)
        self.cStar = self.CEA.get_Cstar(Pc=self.injectionPressure, MR=self.mixtureRatio)
        self.thrustCoefficient = self.CEA.get_PambCf(Pc=self.injectionPressure, MR=self.mixtureRatio, eps=self.expansionRatio, Pamb=self.ambientPressure)[0]
        self.speciesMassFractions = self.CEA.get_SpeciesMassFractions(Pc=self.injectionPressure, MR=self.mixtureRatio, eps=self.expansionRatio, frozen=c, frozenAtThroat=tc, min_fraction=0.005)
        
        temperatures = self.CEA.get_Temperatures(Pc=self.injectionPressure, MR=self.mixtureRatio, eps=self.expansionRatio, frozen=c, frozenAtThroat=tc)
        self.chamberTemp = temperatures[0]
        self.throatTemp = temperatures[1]
        self.exitTemp = temperatures[2]
    

        throatMolWtGamma = self.CEA.get_Throat_MolWt_gamma(Pc=self.injectionPressure, MR=self.mixtureRatio, eps=self.expansionRatio, frozen=c)
        self.throatMolWt = throatMolWtGamma[0]
        self.throatGamma = throatMolWtGamma[1]
        exitMolWtGamma = self.CEA.get_exit_MolWt_gamma(Pc=self.injectionPressure, MR=self.mixtureRatio, eps=self.expansionRatio)
        self.exitMolWt = exitMolWtGamma[0]
        self.exitGamma = exitMolWtGamma[1]
        chamberMolWtGamma = self.CEA.get_Chamber_MolWt_gamma(Pc=self.injectionPressure, MR=self.mixtureRatio, eps=self.expansionRatio)
        self.chamberMolWt = chamberMolWtGamma[0]
        self.chamberGamma = chamberMolWtGamma[1]
        
        chamberTransport = self.CEA.get_Chamber_Transport(Pc=self.injectionPressure, MR=self.mixtureRatio, eps=self.expansionRatio, frozen=c)
        self.chamberHeatCapacity = chamberTransport[0] 
        self.chamberViscosity = chamberTransport[1] * 10
        self.chamberThermalConductivity = chamberTransport[2] * 100
        self.chamberPrandtlNumber = chamberTransport[3]
        exitTransport = self.CEA.get_Exit_Transport(Pc=self.injectionPressure, MR=self.mixtureRatio, eps=self.expansionRatio, frozen=c)
        self.exitHeatCapacity = exitTransport[0] 
        self.exitViscosity = exitTransport[1] * 10
        self.exitThermalConductivity = exitTransport[2] * 100
        self.exitPrandtlNumber = exitTransport[3]

        self.densities = self.CEA.get_Densities(Pc=self.injectionPressure, MR=self.mixtureRatio, eps=self.expansionRatio, frozen=c, frozenAtThroat=tc)
        self.temperatures = self.CEA.get_Temperatures(Pc=self.injectionPressure, MR=self.mixtureRatio, eps=self.expansionRatio, frozen=c, frozenAtThroat=tc)
        self.exitMachNumber = self.CEA.get_MachNumber(Pc=self.injectionPressure, MR=self.mixtureRatio, eps=self.expansionRatio, frozen=c, frozenAtThroat=tc)
    
    def getExitTransportPropertiesAtExpansionRatio(self, expansionRatio):

        c, tc = self.convertConditionToCEAFlag(self.condition, self.throatCondition)

        return self.CEA.get_Exit_Transport(Pc=self.injectionPressure, MR=self.mixtureRatio, eps=expansionRatio, frozen=c)

    # If the mixture ratio override is not set, then it is assumed maximum specific impulse is wanted
    # This uses CEA to calculate the Isp at intervals specified by searchResolution, starting from a mixture ratio specified by startSearch 
    # (can be increased if you know it'll be above 1 to reduce number of calculations)
    def getMaxSpecificImpulseMixtureRatio(self, searchResolution=0.1, startSearch=1):
        
        # Initialises variables for following loop
        i = startSearch
        I_sp = 1
        I_sp_prev = 0
        expansionRatio = 0

        c, tc = self.convertConditionToCEAFlag(self.condition, self.throatCondition)

        # Performs mixture ratio iterations until the specific impulse starts dropping
        while I_sp >= I_sp_prev:

            I_sp_prev = I_sp

            expansionRatio = self.CEA.get_eps_at_PcOvPe(Pc=self.injectionPressure, MR=i, PcOvPe=self.chamberPressure/self.ambientPressure, frozen=c, frozenAtThroat=tc)
            I_sp = self.CEA.estimate_Ambient_Isp(Pc=self.injectionPressure, MR=i, eps=expansionRatio, Pamb=self.ambientPressure, frozen=c, frozenAtThroat=tc)[0]

            i = i + searchResolution
        
        # Rewinds back one iteration to get the ideal specific impulse mixture ratio
        mixtureRatio = i - searchResolution

        return mixtureRatio

    # rocketCEA uses 0 to indicate equilibrium and 1 to indcate frozen
    # This returns the relevant number from the strings of 'equilibrium' and 'frozen'
    def convertConditionToCEAFlag(self, condition, throatCondition):

        # Converts equlibirum/frozen statements into a flag for the CEA functions
        if condition.lower() == 'equilibrium':
            c = 0
        elif condition.lower() == 'frozen':
            c = 1
        else:
            sys.exit('CEA condition can only be "equilibrium" or "frozen"')
        
        if throatCondition.lower() == 'equilibrium':
            ct = 0
        elif throatCondition.lower() == 'frozen':
            ct = 1
        else:
            sys.exit('CEA throat condition can only be "equilibrium" or "frozen"')
        
        return c, ct

    def getMassFlowRate(self):

        self.propellantMassFlowRate = self.thrust / (self.specificImpulse * G)

    def getExitVelocity(self):

        RSpecific = R * 1000 / self.exitMolWt
        self.exitVelocity = sqrt(self.exitGamma * RSpecific * self.temperatures[-1]) * self.exitMachNumber

    def getExitSizes(self):

        self.exitArea = self.propellantMassFlowRate / (self.exitVelocity * self.densities[-1])
        self.exitRadius = areaToRadius(self.exitArea)
        self.exitDiameter = self.exitRadius * 2

    def getThroatSizes(self):

        self.throatArea = self.exitArea / self.expansionRatio
        self.throatRadius = areaToRadius(self.throatArea)
        self.throatDiameter = self.throatRadius * 2

    def getChamberGeometry(self, lStar, contractionLength, entranceRadiusOfCurvatureFactor=1.5, throatEntranceStartAngle=(-135), numberOfPoints=100):

        self.combustionChamber = CombustionChamber(lStar, 
            self.throatRadius, 
            self.CR, 
            contractionLength, 
            entranceRadiusOfCurvatureFactor=entranceRadiusOfCurvatureFactor, 
            throatEntranceStartAngle=throatEntranceStartAngle, 
            numberOfPoints=numberOfPoints)

    def getRaoBellNozzleGeometry(self, lengthFraction, numberOfPoints=300):

        self.nozzle = RaoBellNozzle(self.expansionRatio, self.throatRadius, lengthFraction, numberOfPoints=numberOfPoints)

    def getConicalNozzleGeometry(self, divergenceHalfAngle=15, numberOfPoints=300):

        self.nozzle = ConicalNozzle(self.expansionRatio, self.throatRadius, divergenceHalfAngle=divergenceHalfAngle, numberOfPoints=numberOfPoints)

    def getThrustChamberCoords(self):

        self.throatAverageRadiusOfCurvature = (self.nozzle.throatRadiusOfCurvature + self.combustionChamber.entranceRadiusOfCurvature) / 2
        self.axialCoords = np.concatenate((self.combustionChamber.axialCoords, self.nozzle.axialCoords))
        self.radialCoords = np.concatenate((self.combustionChamber.radialCoords, self.nozzle.radialCoords))

    def plotGeometry(self, part="thrust chamber", show=True, save=False, plotColor="black"):

        fig, ax = plt.subplots()

        if part.lower() == "thrust chamber":

            self.getThrustChamberCoords()
            ax.plot(self.axialCoords, self.radialCoords, color=plotColor)

        elif part.lower() == "combustion chamber":

            ax.plot(self.combustionChamber.axialCoords, self.combustionChamber.radialCoords, color=plotColor)

        elif part.lower() == "nozzle":

            ax.plot(self.nozzle.axialCoords, self.nozzle.radialCoords, color=plotColor)

        else:

            sys.exit("Please give a valid part to plot")

        plt.axis('square')
        plt.ylim(bottom=0)

        if show == True:

            plt.show()
