from PropTools.SubSystems.Engine.Propellant.propellant import Propellant
from PropTools.SubSystems.Engine.Tank.tank import Tank
from PropTools.SubSystems.Engine.Turbopump.pump import Pump
from PropTools.SubSystems.Engine.Turbopump.turbine import Turbine
from PropTools.SubSystems.Engine.ThrustChamber.thrustChamber import ThrustChamber
from PropTools.SubSystems.Engine.ThrustChamber.regenerativeCooling import RegenerativeCooling

class Line:

    def __init__(self, tank: Tank, intermediateComponents, thrustChamber: ThrustChamber):
        
        propellantState = Propellant(tank.outletState.name)


class Cycle:

    def __init__(self, fuelLine, oxLine):

        pass

