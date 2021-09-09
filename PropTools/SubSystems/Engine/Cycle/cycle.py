from copy import deepcopy
from PropTools.SubSystems.Engine.Cycle.component import Component
from PropTools.SubSystems.Engine.Propellant.propellant import Propellant
from PropTools.SubSystems.Engine.Tank.tank import Tank
from PropTools.SubSystems.Engine.Turbopump.pump import Pump
from PropTools.SubSystems.Engine.Turbopump.turbine import Turbine
from PropTools.SubSystems.Engine.ThrustChamber.thrustChamber import ThrustChamber
from PropTools.SubSystems.Engine.ThrustChamber.regenerativeCooling import RegenerativeCooling

class Line:

    def __init__(self, inletState: Propellant, massFlowRate: float, components: "list[Component]", convergenceCriteria=0.1):

        self.inletState = inletState
        self.massFlowRate = massFlowRate
        self.components = components
        self.convergenceCriteria = convergenceCriteria

        self.stateTracker = deepcopy(self.inletState)
        self.states = [deepcopy(self.stateTracker)]

        for i in components:

            if type(i) is Pump:

                i.calculate(self.stateTracker, self.massFlowRate)

            elif type(i) is Turbine:

                i.calculate(self.stateTracker, self.massFlowRate)

            elif type(i) is RegenerativeCooling:

                i.calculate(self.stateTracker, self.massFlowRate, convergenceCriteria=self.convergenceCriteria)


            self.stateTracker = deepcopy(i.outletState)
            self.states.append(deepcopy(self.stateTracker))

        self.outletState = self.states[-1]

class Cycle:

    def __init__(self, fuelLine: Line, oxLine: Line, thrustChamber: ThrustChamber):

        self.fuelLine = fuelLine
        self.oxLine = oxLine
        self.thrustChamber = thrustChamber

        self.powerBalance = 0
        self.pumpPower = 0
        self.turbinePower = 0

        for i in fuelLine.components:

            if type(i) is Pump:

                self.pumpPower += i.power

            elif type(i) is Turbine:

                self.turbinePower += i.power

        for i in oxLine.components:

            if type(i) is Pump:

                self.pumpPower += i.power

            elif type(i) is Turbine:

                self.turbinePower += i.power

        self.powerBalance = self.turbinePower + self.pumpPower

        print("---")
        print("Total Pump Power: " + str(self.pumpPower) + " W")
        print("Total Turbine Power: " + str(self.turbinePower) + " W")
        print("Cycle Power Balance: " + str(self.powerBalance) + " W")

        

            




        



        


