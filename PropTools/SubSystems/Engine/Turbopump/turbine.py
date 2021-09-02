from PropTools.SubSystems.Engine.Cycle.component import Component
from PropTools.SubSystems.Engine.Propellant.propellant import Propellant

class Turbine(Component):

    def __init__(self, inletState: Propellant, outletPressure, massFlowRate, isentropicEfficiency):

        super().__init__()
        self.inletState = inletState
        self.outletPressure = outletPressure
        self.isentropicEfficiency = isentropicEfficiency

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
