import csv
import sys
from math import log
import LRETools.SubSystems.Body.airframe as airframe
import LRETools.SubSystems.Engine.engine as engine
import LRETools.SubSystems.FlightSimulation.flight as flight
import LRETools.SubSystems.Engine.Cycle.propellant as propellant

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

    def flightSimulation(self, timeStep=0.1, variablesToPlot="All", showPlot=False, saveFileName="Simulation.png", plotAbsoluteValues=True):
        
        self.flightSimulation = flight.Flight(self.airframe.frontalArea, self.airframe.cdValues, self.airframe.machValues)
        self.flightSimulation.flightSimulation(self.wetMass, self.propellantMass, self.engine.thrust, self.engine.specificImpulse, timeStep=timeStep, saveAbsoluteValues=plotAbsoluteValues)
        self.flightSimulation.plotGraph(variablesToPlot=variablesToPlot, show=showPlot, fileName=saveFileName)