from PropTools.SubSystems.Engine.Tank.tank import Tank
from PropTools.SubSystems.Engine.Turbopump.pump import Pump
from PropTools.SubSystems.Engine.Turbopump.turbine import Turbine
from PropTools.SubSystems.Engine.ThrustChamber.thrustChamber import ThrustChamber
from PropTools.SubSystems.Engine.ThrustChamber.regenerativeCooling import RegenerativeCooling

class Line:

    def __init__(self, tank, intermediateComponents, thrustChamber):

        for i in intermediateComponents:

            if i.lower() == "pump":

                pass

            elif i.lower() == "turbine":

                pass

            elif i.lower() == "regenerative cooling" or "regen cooling" or "regen" or "hx":

                pass

class Cycle:

    def __init__(self, fuelLine, oxLine):

        pass

