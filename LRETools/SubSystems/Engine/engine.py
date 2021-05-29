import sys
from rocketcea.cea_obj_w_units import CEA_Obj
from LRETools.SubSystems.Engine.Cycle.propellant import Propellant
from LRETools.SubSystems.Engine.Cycle.cycle import *

class Engine:

    def __init__(self, fuelName, oxName):
        
        self.fuel = Propellant(fuelName)
        self.ox = Propellant(oxName)

        self.cycle = ExpanderCycle(self.fuel.name, self.ox.name)

        print(self.fuel.ceaName)
        print(self.cycle.fuel.ceaName)

        # Initialise variables that will be defined using other functions
        self.mass = None
        self.chamberPressure = None
        self.thrust = None
        self.expansionRatio = None
        self.specificImpulse = None
        self.chamberTemp = None
        self.throatMolWtGamma = None
        self.throatMolWt = None
        self.gamma = None
        self.propellantMassFlowRate = None
        self.fuelMassFlowRate = None
        self.oxMassFlowRate = None

    def setMass(self, mass):

        self.mass = mass

    def setChamberPressure(self, chamberPressure):

        self.chamberPressure = chamberPressure

    def setThrust(self, thrust):

        self.thrust = thrust

    # Iterates through mixture ratios, performing CEA on each one, to find the ratio that produces the highest specific impulse
    def findIdealMixtureRatio(self, ambientPressure=1.01325, condition='equilibrium', throatCondition='equilibrium', searchResolution=0.1, startSearch=2):

        CEA = CEA_Obj(oxName=self.ox.ceaName, fuelName=self.fuel.ceaName)

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

            expansionRatio = CEA.get_eps_at_PcOvPe(Pc=self.chamberPressure, MR=i, PcOvPe=self.chamberPressure/ambientPressure, frozen=c, frozenAtThroat=ct)
            I_sp = CEA.get_Isp(Pc=self.chamberPressure, MR=i, eps=expansionRatio, frozen=c, frozenAtThroat=ct)

            i = i + searchResolution
        
        # Rewinds back one iteration to get the ideal specific impulse mixture ratio
        self.idealMixtureRatio = i - searchResolution

        # Gets useful parameters at the ideal mixture ratio
        self.expansionRatio = CEA.get_eps_at_PcOvPe(Pc=self.chamberPressure, MR=self.idealMixtureRatio, PcOvPe=self.chamberPressure/ambientPressure, frozen=c, frozenAtThroat=ct)
        self.specificImpulse = CEA.get_Isp(Pc=self.chamberPressure, MR=self.idealMixtureRatio, eps=expansionRatio, frozen=c, frozenAtThroat=ct)
        self.chamberTemp = CEA.get_Tcomb(Pc=self.chamberPressure, MR=i)
        self.speciesMassFractions = CEA.get_SpeciesMassFractions(Pc=self.chamberPressure, MR=self.idealMixtureRatio, eps=expansionRatio, frozen=c, frozenAtThroat=ct, min_fraction=0.005)
        throatMolWtGamma = CEA.get_Throat_MolWt_gamma(Pc=self.chamberPressure, MR=self.idealMixtureRatio, eps=expansionRatio, frozen=ct)
        self.throatMolWt = throatMolWtGamma[0]
        self.gamma = throatMolWtGamma[1]

        self.cycle.setChamberConditions(self.chamberPressure, self.chamberTemp)

    def findPropellantMassFlowRates(self):

        self.propellantMassFlowRate = self.thrust / (self.specificImpulse * 9.81)
        self.fuelMassFlowRate = self.propellantMassFlowRate / self.idealMixtureRatio
        self.oxMassFlowRate = self.fuelMassFlowRate * self.idealMixtureRatio

        self.cycle.setMassFlowRate(self.fuelMassFlowRate, self.oxMassFlowRate)


    