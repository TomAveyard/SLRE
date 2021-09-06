from PropTools.SubSystems.Engine.Cycle.component import Component
from PropTools.SubSystems.Engine.Propellant.propellant import Propellant

class Turbine(Component):

    def __init__(self, isentropicEfficiency, pressureLoss=None, outletPressure=None):

        super().__init__() # Contains self.inletState, self.outletState, self.type
        self.isentropicEfficiency = isentropicEfficiency
        self.type = "turbine"

        if pressureLoss == None and outletPressure == None:
            exit("Please specify either a pressure loss or an outlet pressure for the turbine")
        elif pressureLoss != None and outletPressure != None:
            exit("Please specify either a pressure loss or an outlet pressure for the turbine, not both")

        self.pressureLoss = pressureLoss
        self.outletPressure = outletPressure

        self.massFlowRate = None
        self.deltaHIsentropic = None
        self.deltaHReal = None
        self.massFlowRate = None
        self.power = None

    def calculate(self, inletState: Propellant, massFlowRate):

        self.inletState = inletState

        if self.pressureLoss != None:
            self.outletPressure = self.inletState.P - self.pressureLoss

        # Defines the isentropic outlet state
        outletStateIsentropic = Propellant(self.inletState.name)
        outletStateIsentropic.defineState("S", self.inletState.S, "P", self.outletPressure)

        # Calculates the isentropic enthalpy change
        self.deltaHIsentropic = self.inletState.H - outletStateIsentropic.H

        # Calculates the real enthalpy change
        self.deltaHReal = self.deltaHIsentropic * self.isentropicEfficiency

        # Defines real state at the outlet
        self.outletState = Propellant(self.inletState.name)
        self.outletState.defineState("H", inletState.H - self.deltaHReal, "P", self.outletPressure)

        self.massFlowRate = massFlowRate
        self.power = self.deltaHReal * self.massFlowRate
