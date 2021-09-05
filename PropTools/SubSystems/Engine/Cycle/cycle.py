from re import L
from PropTools.SubSystems.Engine.Cycle.component import Component
from PropTools.SubSystems.Engine.Propellant.propellant import Propellant
from PropTools.SubSystems.Engine.Tank.tank import Tank
from PropTools.SubSystems.Engine.Turbopump.pump import Pump
from PropTools.SubSystems.Engine.Turbopump.turbine import Turbine
from PropTools.SubSystems.Engine.ThrustChamber.thrustChamber import ThrustChamber
from PropTools.SubSystems.Engine.ThrustChamber.regenerativeCooling import RegenerativeCooling

class Line:

    def __init__(self, tank: Tank, intermediateComponents: "list[Component]", thrustChamber: ThrustChamber, propellant):
        
        propellantState = Propellant(tank.outletState.name)
        propellantState.defineState("T", tank.outletState.T, "P", tank.outletState.P)

        self.states = [propellantState]

        if propellant.lower() == "fuel":
            self.massFlowRate = thrustChamber.fuelMassFlowRate
        elif propellant.lower() == "ox" or propellant.lower() == "oxidiser" or propellant.lower() == "oxidizer":
            self.massFlowRate = thrustChamber.oxMassFlowRate

        if "regenerative cooling" in intermediateComponents:

            for i in intermediateComponents:

                if i.type == "pump":

                    pass

                if i.type == "turbine":

                    pass

                if i.type == "regenerative cooling":

                    pass

        else:

            for i in intermediateComponents:

                if i.type == "pump":

                    i.calculate(propellantState, thrustChamber.injectionPressure, self.massFlowRate)
                    self.states.append(i.outletState)
                    propellantState.defineState("T", i.outletState.T, "P", i.outletState.P)

                if i.type == "turbine":

                    pass

        

class Cycle:

    def __init__(self, fuelLine, oxLine):

        pass
