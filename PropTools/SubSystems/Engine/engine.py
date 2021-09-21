import sys
from rocketcea.cea_obj_w_units import CEA_Obj
from PropTools.SubSystems.Engine.Propellant.propellant import Propellant
from PropTools.SubSystems.Engine.Cycle.cycle import *

#TODO: Update to reflect changes to underlying architecture

class Engine:

    def __init__(self, fuelName, oxName):
        
        self.fuel = Propellant(fuelName)
        self.ox = Propellant(oxName)

        #self.cycle = ExpanderCycle(self.fuel.name, self.ox.name)

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

    def findPropellantMassFlowRates(self):

        self.propellantMassFlowRate = self.thrust / (self.specificImpulse * 9.81)
        self.fuelMassFlowRate = self.propellantMassFlowRate / self.idealMixtureRatio
        self.oxMassFlowRate = self.fuelMassFlowRate * self.idealMixtureRatio

        self.cycle.setMassFlowRate(self.fuelMassFlowRate, self.oxMassFlowRate)


    