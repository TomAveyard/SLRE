from PropTools.SubSystems.Engine.Propellant.propellant import Propellant
from PropTools.SubSystems.Engine.ThrustChamber import thrustChamber, combustionChamber, nozzle, regenerativeCooling

#TODO: Update to reflect changes to underlying architecture

class Engine:

    def __init__(self, fuelName, oxName):
        
        self.fuel = Propellant(fuelName)
        self.ox = Propellant(oxName)

    


    