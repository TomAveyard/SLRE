from PropTools.SubSystems.Engine.Propellant.propellant import Propellant
from PropTools.SubSystems.Engine.ThrustChamber.thrustChamber import ThrustChamber
from PropTools.SubSystems.Engine.Cycle.cycle import Cycle

#Wrapper class that uses underlying objects to represent the engine as a whole

class Engine:

    def __init__(self, thrustChamber: ThrustChamber = None, cycle: Cycle = None, mass: float = None):

        self.thrustChamber = thrustChamber
        self.cycle = cycle
        self.mass = mass


        
        

    




    


    