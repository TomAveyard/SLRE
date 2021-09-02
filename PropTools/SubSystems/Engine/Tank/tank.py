from PropTools.SubSystems.Engine.Cycle.component import Component
from PropTools.SubSystems.Engine.Propellant.propellant import Propellant

# Defines state of propellant in the tank from the stated temperature and pressure of the tank
class Tank:

    def __init__(self, outletState:Propellant):

        self.outletState = outletState