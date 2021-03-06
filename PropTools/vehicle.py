import csv
import sys
from math import log
import PropTools.SubSystems.Body.airframe as airframe
import PropTools.SubSystems.Engine.engine as engine

class Vehicle:

    def __init__(self, name, fuel, ox):

        self.name = name
        self.G = 9.81
        self.airframe = airframe.Airframe()
        self.engine = engine.Engine(fuel, ox)

        # Initialise variables that will be defined using other functions
        self.fuelMass = None
        self.oxMass = None
        self.propellantMass = None
        self.dryMass = None
        self.wetMass = None
        self.deltaV = None
        self.flightSimulation = None
        
    def calculateMasses(self):

        self.fuelMass = self.airframe.totalTankVolume / ((self.engine.idealMixtureRatio / self.engine.cycle.oxTank.D) + (1 / self.engine.cycle.fuelTank.D))
        self.oxMass = self.fuelMass * self.engine.idealMixtureRatio
        self.propellantMass = self.fuelMass + self.oxMass
        self.dryMass = self.airframe.mass + self.engine.mass
        self.wetMass = self.dryMass + self.propellantMass  

    def calculateDeltaV(self):
        
        self.deltaV = self.engine.specificImpulse * self.G * log(self.wetMass / self.dryMass)