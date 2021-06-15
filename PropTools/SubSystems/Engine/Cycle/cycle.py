import sys
from PropTools.SubSystems.Engine.Cycle.propellant import Propellant
from PropTools.SubSystems.Engine.Cycle.cyclecomponents import *

class ExpanderCycle:

    def __init__(self, fuel, ox):
        
        self.fuel = Propellant(fuel)
        self.ox = Propellant(ox)

        # Initialise variables that will be defined using other functions
        self.fuelTank = None
        self.oxTank = None
        self.chamberPressure = None
        self.chamberTemperature = None
        self.injectorPressureLoss = None
        self.fuelMassFlowRate = None
        self.oxMassFlowRate = None
        self.fuelPumpEfficiency = None
        self.oxPumpEfficiency = None
        self.turbineEfficiency = None
        self.HXHeatTransferRate = None
        self.HXPressureLoss = None
        self.T_g = None
        self.T_l = None
        self.h_g = None
        self.h_l = None
        self.t_w = None
        self.k = None
        self.A = None
        self.oxPump = None
        self.fuelPump = None
        self.injectorInletPressure = None
        self.regenCooling = None
        self.turbine = None

    
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
    def setRegenCoolingEstimate(self, h_g, h_l, T_g, T_l, t_w, k, A, pressureLoss):
        
        # h_g = gas-side coefficient
        self.h_g = h_g
        # h_l = liquid-side coefficient
        self.h_l = h_l
        # T_g = gas-side temperature
        self.T_g = T_g
        # T_l = liquid-side temperature. The value provided here is an initial guess
        self.T_l = T_l
        # t_w = wall thickness
        self.t_w = t_w
        # k = conductivity
        self.k = k
        # A = area over which heat transfer takes place
        self.A = A

        heatTransferRatePerArea = (self.T_g - self.T_l) / ((1 / self.h_g) + (self.t_w / self.k) + (1 / self.h_l))
        heatTransferRate = heatTransferRatePerArea * self.A

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

                averageRegenCoolingTemp = (self.regenCooling.T + self.fuelPump.T) / 2
                self.T_l = averageRegenCoolingTemp
                # Recalculates rate of heat transfer
                self.setRegenCoolingEstimate(self.h_g, self.h_l, self.T_g, self.T_l, self.t_w, self.k, self.A, self.HXPressureLoss)

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
