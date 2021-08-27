from PropTools.SubSystems.Engine.Propellant.propellant import Propellant

# Calculates state of propellant at a turbine outlet from an inlet state, target outlet pressure, and isentropic efficiency
class Turbine(Propellant):

    def __init__(self, propellant, inletState, outletPressure, isentropicEfficiency):
        
        # Initialises super class
        super().__init__(propellant)

        self.isentropicEfficiency = isentropicEfficiency

        # Defines the isentropic outlet state
        self.defineState("S", inletState.S, "P", outletPressure)
        # Calculates the isentropic enthalpy change
        self.deltaHIsentropic = inletState.H - self.H

        # Defines real state at the outlet
        self.deltaHReal = self.deltaHIsentropic * self.isentropicEfficiency

        # Defines real state at the outlet
        self.defineState("H", inletState.H - self.deltaHReal, "P", outletPressure)

        # Initialise variables that will be defined using other functions
        self.massFlowRate = None
        self.power = None

    # Calculates the power required by the pump for a given mass flow rate
    def calculatePower(self, massFlowRate):

        self.massFlowRate = massFlowRate
        self.power = self.deltaHReal * massFlowRate