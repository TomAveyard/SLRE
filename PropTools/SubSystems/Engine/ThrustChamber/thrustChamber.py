import sys
from rocketcea.cea_obj_w_units import CEA_Obj
from PropTools.SubSystems.Engine.Cycle.Propellant.propellant import Propellant
from PropTools.SubSystems.Engine.ThrustChamber.combustionChamber import CombustionChamber
from PropTools.SubSystems.Engine.ThrustChamber.nozzle import Nozzle
import matplotlib.pyplot as plt

class ThrustChamber:

    def __init__(self, 
                fuelName, 
                oxName, 
                chamberPressure, 
                mixtureRatioOverride=False, 
                ambientPressure=1.01325, 
                facCR=None,
                facPlenumPressureSpecified=True,
                mixtureRatioSearchResolution=0.1, 
                mixtureRatioSearchStart=1, 
                condition='equilibrium', 
                throatCondition='equilibrium'):

        self.fuel = Propellant(fuelName)
        self.ox = Propellant(oxName)

        self.chamberPressure = chamberPressure
        self.ambientPressure = ambientPressure

        self.facCR = facCR

        # If a finite area combustor contraction ratio is specified and the chamber pressure given is the plenum/combustion pressure
        # then this estimates an injector face pressure for use in the calculations.
        # Otherwise the injection pressure is the chamber pressure specified.
        if self.facCR != None and facPlenumPressureSpecified == True:
            self.facPlenumPressureSpecified = facPlenumPressureSpecified
            self.injectionPressureOverCombustionPressure = 1.0 + 0.54 / self.facCR**2.2
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

    # Returns a cea object using the unit system seen below
    def getCEAObject(self):

        # If facCR (finite area combustor contraction ratio) is not specified, then an infinite area combustor is assumed
        if self.facCR == None:

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

        else:

            CEAObject = CEA_Obj(oxName=self.ox.ceaName, 
                fuelName=self.fuel.ceaName,
                fac_CR=self.facCR, 
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
        self.chamberTemp = self.CEA.get_Tcomb(Pc=self.injectionPressure, MR=self.mixtureRatio)
        self.speciesMassFractions = self.CEA.get_SpeciesMassFractions(Pc=self.injectionPressure, MR=self.mixtureRatio, eps=self.expansionRatio, frozen=c, frozenAtThroat=tc, min_fraction=0.005)
        throatMolWtGamma = self.CEA.get_Throat_MolWt_gamma(Pc=self.injectionPressure, MR=self.mixtureRatio, eps=self.expansionRatio, frozen=tc)
        self.throatMolWt = throatMolWtGamma[0]
        self.gamma = throatMolWtGamma[1]

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
