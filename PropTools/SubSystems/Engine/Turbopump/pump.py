from sys import exit
from PropTools.SubSystems.Engine.Cycle.component import Component
from PropTools.SubSystems.Engine.Propellant.propellant import Propellant

class Pump(Component):

    def __init__(self, isentropicEfficiency, pressureRise=None, outletPressure=None):

        super().__init__() # Contains self.inletState, self.outletState, self.type
        self.isentropicEfficiency = isentropicEfficiency
        self.type = "pump"

        if pressureRise == None and outletPressure == None:
            exit("Please specify either a pressure rise or an outlet pressure for the pump")
        elif pressureRise != None and outletPressure != None:
            exit("Please specify either a pressure rise or an outlet pressure for the pump, not both")

        self.pressureRise = pressureRise
        self.outletPressure = outletPressure

        self.massFlowRate = None
        self.deltaHIsentropic = None
        self.deltaHReal = None
        self.massFlowRate = None
        self.power = None

    def calculate(self, inletState: Propellant, massFlowRate):

        self.inletState = inletState

        if self.pressureRise != None:
            self.outletPressure = self.inletState.P + self.pressureRise

        # Defines the isentropic outlet state
        outletStateIsentropic = Propellant(self.inletState.name)
        outletStateIsentropic.defineState("S", self.inletState.S, "P", self.outletPressure)

        # Calculates the isentropic enthalpy change
        self.deltaHIsentropic = self.inletState.H - outletStateIsentropic.H

        # Calculates the real enthalpy change
        self.deltaHReal = self.deltaHIsentropic / self.isentropicEfficiency

        # Defines real state at the outlet
        self.outletState = Propellant(self.inletState.name)
        self.outletState.defineState("H", inletState.H - self.deltaHReal, "P", self.outletPressure)

        self.massFlowRate = massFlowRate
        self.power = self.deltaHReal * self.massFlowRate
