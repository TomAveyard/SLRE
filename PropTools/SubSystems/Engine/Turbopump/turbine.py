from PropTools.SubSystems.Engine.Cycle.component import Component
from PropTools.SubSystems.Engine.Propellant.propellant import Propellant

class Turbine(Component):

    def __init__(self, inletState: Propellant, outletPressure, massFlowRate, isentropicEfficiency):

        super().__init__() # Contains self.inletState, self.outletState, self.type
        self.isentropicEfficiency = isentropicEfficiency
        self.type = "Turbine"

        self.deltaHIsentropic = None
        self.deltaHReal = None
        self.massFlowRate = None
        self.power = None

    def calculate(self, inletState: Propellant, outletPressure, massFlowRate):

        self.inletState = inletState
        outletPressure = outletPressure

        # Defines the isentropic outlet state
        outletStateIsentropic = Propellant(self.inletState.name)
        outletStateIsentropic.defineState("S", self.inletState.S, "P", outletPressure)

        # Calculates the isentropic enthalpy change
        self.deltaHIsentropic = self.inletState.H - outletStateIsentropic.H

        # Calculates the real enthalpy change
        self.deltaHReal = self.deltaHIsentropic * self.isentropicEfficiency

        # Defines real state at the outlet
        self.outletState = Propellant(self.inletState.name)
        self.outletState.defineState("H", inletState.H - self.deltaHReal, "P", outletPressure)

        # Initialise variables that will be defined using other functions
        self.massFlowRate = massFlowRate
        self.power = self.deltaHReal * massFlowRate
