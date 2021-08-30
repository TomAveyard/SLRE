from PropTools.SubSystems.Engine.Propellant.propellant import Propellant

# Defines state of propellant in the tank from the stated temperature and pressure of the tank
class Tank(Propellant):

    def __init__(self, propellant, temp, pressure):
        
        # Initialises super class
        super().__init__(propellant)
        self.defineState("T", temp, "P", pressure)

        self.totalTankVolume = None