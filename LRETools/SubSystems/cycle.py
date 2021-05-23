import sys
from LRETools.SubSystems.propellant import Propellant
from CoolProp.CoolProp import PropsSI

# Defines state of propellant in the tank from the stated temperature and pressure of the tank
class Tank(Propellant):

    def __init__(self, propellant, temp, pressure):
        
        # Initialises super class
        super().__init__(propellant)
        self.defineState("T", temp, "P", pressure)

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

class ExpanderCycle:

    def __init__(self, fuel, ox):
        
        self.fuel = Propellant(fuel)
        self.ox = Propellant(ox)
    
    def setFuelTank(self, fuelTemp, fuelPressure):

        self.fuelTank = Tank(self.fuel.name, fuelTemp, fuelPressure)

    def setOxTank(self, oxTemp, oxPressure):

        self.oxTank = Tank(self.ox.name, oxTemp, oxPressure)

    def setChamberConditions(self, chamberPressure, chamberTemperature):

        self.chamberPressure = chamberPressure
        self.chamberTemperature = chamberTemperature

    def setInjector(self, injectorPressureLoss):

        self.injectorPressureLoss = injectorPressureLoss

    def setMassFlowRate(self, fuelMassFlowRate, oxMassFlowRate):

        self.fuelMassFlowRate = fuelMassFlowRate
        self.oxMassFlowRate = oxMassFlowRate

    def setFuelPumpEfficiency(self, efficiency):

        self.fuelPumpEfficiency = efficiency

    def setOxPumpEfficiency(self, efficiency):

        self.oxPumpEfficiency = efficiency

    def setTurbineEfficiency(self, efficiency):

        self.turbineEfficiency = efficiency

    # If the rate of heat transfer is known, then it can be set as a constant
    def setRegenCoolingConstant(self, heatTransferRateConstant, pressureLoss):

        self.HXHeatTransferRate = heatTransferRateConstant
        self.HXPressureLoss = pressureLoss

    # If the rate of heat transfer is not known, it can be estimated with 7 parameters set in this function
    def setRegenCoolingEstimate(self, T_g, T_l, h_g, t_w, k, h_l, A, pressureLoss):

        # T_g = gas-side temperature
        self.T_g = T_g
        # T_l = liquid-side temperature. The value provided here is an initial guess
        self.T_l = T_l
        self.h_g = h_g
        self.t_w = t_w
        self.k = k
        self.h_l = h_l
        # A = area over which heat transfer takes place
        self.A = A

        heatTransferRatePerArea = (T_g - T_l) / ((1 / h_g) + (t_w / k) + (1 / h_l))
        heatTransferRate = heatTransferRatePerArea * A

        self.HXHeatTransferRate = heatTransferRate
        self.HXPressureLoss = pressureLoss

    def solveCycle(self, searchInterval=100000, accuracy=1, searchDivision=2, updateHeatTransferEstimate=True):

        # Counter for number of iterations
        count = 0

        print("Solving...")

        # Sets the pressure required at the injector inlet
        self.injectorInletPressure = self.chamberPressure + self.injectorPressureLoss

        # Calculates oxidiser pump outlet state and power
        self.oxPump = Pump(self.ox.name, self.oxTank, self.injectorInletPressure, self.oxPumpEfficiency)
        self.oxPump.calculatePower(self.oxMassFlowRate)

        # Initialises turbine power variable used for the loop exit condition
        turbinePower = 0

        # Estimates an outlet pressure for the fuel pump to get initialise the iterative loop
        fuelPumpOutletPressureGuess = self.injectorInletPressure + self.HXPressureLoss + 1e5

        # Estimation for total pump power to initialise the iterative loop
        totalPumpPower = self.oxPump.power + ((fuelPumpOutletPressureGuess * self.fuelMassFlowRate) / (self.fuelTank.D * self.fuelPumpEfficiency))

        # Loop will exit when the power required for both pumps is balanced with the power provided by the turbine, to a specified accuracy
        while round(totalPumpPower, -accuracy) != round(-turbinePower, -accuracy):
            
            # Increases the guessed fuel pump outlet pressure by a specified amount each iteration
            fuelPumpOutletPressureGuess = fuelPumpOutletPressureGuess + searchInterval

            # Calculates fuel pump outlet state and power
            self.fuelPump = Pump(self.fuel.name, self.fuelTank, fuelPumpOutletPressureGuess, self.fuelPumpEfficiency)
            self.fuelPump.calculatePower(self.fuelMassFlowRate)

            totalPumpPower = self.fuelPump.power + self.oxPump.power
            
            # Calculates the regenerative cooling loop outlet state
            self.regenCooling = HeatExchanger(self.fuel.name, self.fuelPump, self.HXHeatTransferRate, self.fuelMassFlowRate, self.HXPressureLoss)

            # If we want to update the rate of heat transfer with the calculated average liquid side temperature, the following will do that
            if updateHeatTransferEstimate == True:

                # Checks to see if we are using an estimated heat transfer
                try:
                    self.T_l
                    averageRegenCoolingTemp = (self.regenCooling.T + self.fuelPump.T) / 2
                    self.T_l = averageRegenCoolingTemp
                    # Recalculates rate of heat transfer
                    self.setRegenCoolingEstimate(self.T_g, self.T_l, self.h_g, self.t_w, self.k, self.h_l, self.A, self.HXPressureLoss)
                except:
                    sys.exit("Solution aborted: updateHeatTransferEstimate=True but the heat transfer rate was set as constant")

            # Calculates turbine outlet state and power
            self.turbine = Turbine(self.fuel.name, self.regenCooling, self.injectorInletPressure, self.turbineEfficiency)
            self.turbine.calculatePower(self.fuelMassFlowRate)

            # Updates the exit condition variable
            turbinePower = self.turbine.power
            
            # Adjust searching parameter to achieve higher accuracy
            if -totalPumpPower > self.turbine.power and searchInterval < 0:
                searchInterval = -searchInterval / searchDivision
            elif -totalPumpPower < self.turbine.power and searchInterval > 0:
                searchInterval = -searchInterval / searchDivision

            count = count + 1

            print("Current Total Pump Power: " + str(totalPumpPower) + " | Current Turbine Power: " + str(self.turbine.power), end="\r")

        print("\nSolution Completed in " + str(count) + " iterations")

            









