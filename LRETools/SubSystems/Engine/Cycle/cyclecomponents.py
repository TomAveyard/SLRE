from LRETools.SubSystems.Engine.propellant import Propellant

# Defines state of propellant in the tank from the stated temperature and pressure of the tank
class Tank(Propellant):

    def __init__(self, propellant, temp, pressure):
        
        # Initialises super class
        super().__init__(propellant)
        self.defineState("T", temp, "P", pressure)

        self.totalTankVolume = None

# Calculates state of propellant at a pump outlet from an inlet state, target outlet pressure, and isentropic efficiency
class Pump(Propellant):

    def __init__(self, propellant, inletState, outletPressure, isentropicEfficiency):

        # Initialises super class
        super().__init__(propellant)

        self.isentropicEfficiency = isentropicEfficiency

        # Defines the isentropic outlet state
        self.defineState("S", inletState.S, "P", outletPressure)

        # Calculates the isentropic enthalpy change
        self.deltaHIsentropic = inletState.H - self.H

        # Calculates the real enthalpy change
        self.deltaHReal = self.deltaHIsentropic / self.isentropicEfficiency

        # Defines real state at the outlet
        self.defineState("H", inletState.H - self.deltaHReal, "P", outletPressure)

        # Initialise variables that will be defined using other functions
        self.massFlowRate = None
        self.power = None

    # Calculates the power required by the pump for a given mass flow rate
    def calculatePower(self, massFlowRate):

        self.massFlowRate = massFlowRate
        self.power = self.deltaHReal * massFlowRate

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

class HeatExchanger(Propellant):

    def __init__(self, propellant, inletState, heatTransferRate, massFlowRate, pressureLoss):

        # Initialises super class
        super().__init__(propellant)

        self.heatTransferRate = heatTransferRate
        self.massFlowRate = massFlowRate
        self.pressureLoss = pressureLoss

        # Calculates the enthalpy change from the stated rate of heat transfer
        self.deltaH = self.heatTransferRate / self.massFlowRate

        # Defines state at the outlet of the heat exchanger
        self.defineState("H", inletState.H + self.deltaH, "P", inletState.P - pressureLoss)