import sys
from rocketcea.cea_obj_w_units import CEA_Obj
from PropTools.SubSystems.Engine.Cycle.Propellant.propellant import Propellant

class thrustChamber:

    def __init__(self, fuelName, oxName):

        self.fuel = Propellant(fuelName)
        self.ox = Propellant(oxName)
    
    def setChamberPressure(self, chamberPressure):

        self.chamberPressure = chamberPressure

    def setMixtureRatio(self, mixtureRatio):

        self.mixtureRatio = mixtureRatio

    def getCEA(self, ambientPressure=1.01325, condition='equilibrium', throatCondition='equilibrium'):

        self.CEA = CEA_Obj(oxName=self.ox.ceaName, 
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
            thermal_cond_units = 'W/cm-degC')

    # Iterates through mixture ratios, performing CEA on each one, to find the ratio that produces the highest specific impulse
    def findMaxSpecificImpulseMixtureRatio(self, ambientPressure=1.01325, condition='equilibrium', throatCondition='equilibrium', searchResolution=0.1, startSearch=2):

        self.getCEA()

        # Converts equlibirum/frozen statements into a flag for the CEA functions
        if condition == 'equilibrium':

            c = 0

        elif condition == 'frozen':

            c = 1

        else:

            sys.exit('CEA condition can only be "equilibrium" or "frozen"')
        
        if throatCondition == 'equilibrium':

            ct = 0

        elif throatCondition == 'frozen':

            ct = 1

        else:
            
            sys.exit('CEA throat condition can only be "equilibrium" or "frozen"')

        # Initialises variables for following loop
        i = startSearch
        I_sp = 1
        I_sp_prev = 0
        expansionRatio = 0

        # Performs mixture ratio iterations until the specific impulse starts dropping
        while I_sp >= I_sp_prev:

            I_sp_prev = I_sp

            expansionRatio = self.CEA.get_eps_at_PcOvPe(Pc=self.chamberPressure, MR=i, PcOvPe=self.chamberPressure/ambientPressure, frozen=c, frozenAtThroat=ct)
            I_sp = self.CEA.get_Isp(Pc=self.chamberPressure, MR=i, eps=expansionRatio, frozen=c, frozenAtThroat=ct)

            i = i + searchResolution
        
        # Rewinds back one iteration to get the ideal specific impulse mixture ratio
        self.idealMixtureRatio = i - searchResolution

        # Gets useful parameters at the ideal mixture ratio
        self.expansionRatio = self.CEA.get_eps_at_PcOvPe(Pc=self.chamberPressure, MR=self.idealMixtureRatio, PcOvPe=self.chamberPressure/ambientPressure, frozen=c, frozenAtThroat=ct)
        self.specificImpulse = self.CEA.get_Isp(Pc=self.chamberPressure, MR=self.idealMixtureRatio, eps=expansionRatio, frozen=c, frozenAtThroat=ct)
        self.chamberTemp = self.CEA.get_Tcomb(Pc=self.chamberPressure, MR=i)
        self.speciesMassFractions = self.CEA.get_SpeciesMassFractions(Pc=self.chamberPressure, MR=self.idealMixtureRatio, eps=expansionRatio, frozen=c, frozenAtThroat=ct, min_fraction=0.005)
        throatMolWtGamma = self.CEA.get_Throat_MolWt_gamma(Pc=self.chamberPressure, MR=self.idealMixtureRatio, eps=expansionRatio, frozen=ct)
        self.throatMolWt = throatMolWtGamma[0]
        self.gamma = throatMolWtGamma[1]

    def bartzEquation():
        pass
