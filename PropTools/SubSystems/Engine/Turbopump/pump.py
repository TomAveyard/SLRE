from sys import exit

from fluids.constants.constants import G
from PropTools.SubSystems.Engine.Cycle.component import Component
from PropTools.SubSystems.Engine.Propellant.propellant import Propellant

class Pump(Component):

    def __init__(self, isentropicEfficiency: float = None, pressureRise: float = None, outletPressure: float = None):

        super().__init__() # Contains self.inletState, self.outletState, self.type
        self.isentropicEfficiency = isentropicEfficiency
        self.type = "pump"

        if pressureRise == None and outletPressure == None:
            exit("Please specify either a pressure rise or an outlet pressure for the pump")
        elif pressureRise != None and outletPressure != None:
            exit("Please specify either a pressure rise or an outlet pressure for the pump, not both")

        self.pressureRise = pressureRise
        self.headRise = self.pressureRise / (G * self.inletState.D)
        self.outletPressure = outletPressure

        self.massFlowRate = None
        self.deltaHIsentropic = None
        self.deltaHReal = None
        self.massFlowRate = None
        self.power = None

    def calculate(self, inletState: Propellant = None, massFlowRate: float = None) -> None:

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

    def printResults(self, label="Pump", decimalPlaces: int = 2):

        print(self.printSeperator)
        print(label + " Results")
        print(self.printSeperator)
        print("Substance: " + self.inletState.name)
        print("Mass Flow Rate: " + str(round(self.massFlowRate, decimalPlaces)) + " kg/s")
        print("Inlet Temperature: " + str(round(self.inletState.T, decimalPlaces)) + " K")
        print("Inlet Pressure: " + str(round(self.inletState.P/1e5, decimalPlaces)) + " Bar")
        print("Outlet Temperature: " + str(round(self.outletState.T, decimalPlaces)) + " K")
        print("Outlet Pressure: " + str(round(self.outletState.P/1e5, decimalPlaces)) + " Bar")
        print("Pressure Rise: " + str(round(self.pressureRise, decimalPlaces)) + " Bar")
        print("Head Rise: " + str(round(self.headRise, decimalPlaces)) + " m")
        print("Power: " + str(round(self.power/1e3, decimalPlaces)) + " kW")
        print(self.printSeperator)
