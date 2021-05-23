import csv
import sys
from math import log
import LRETools.SubSystems.cea as cea
import LRETools.SubSystems.flight as flight
import LRETools.SubSystems.propellant as propellant

class Vehicle:

    def __init__(self, name):

        self.name = name
        self.G = 9.81

    def setMasses(self, airFrameMass, engineMass, fuelMass, oxMass):
        
        self.airFrameMass = airFrameMass
        self.engineMass = engineMass
        self.dryMass = airFrameMass + engineMass
        self.fuelMass = fuelMass
        self.oxMass = oxMass
        self.propellantMass = fuelMass + oxMass
        self.wetMass = self.dryMass + self.propellantMass

    def readDragCurve(self, dragCurveCSV):

        machList = []
        cdList = []

        with open(dragCurveCSV) as file:
            reader = csv.reader(file)
            for row in reader:
                machList.append(row[0])
                cdList.append(row[1])

        
        return [machList, cdList]

    def setDragCharacteristics(self, frontalArea, dragCurve):

        self.frontalArea = frontalArea

        if type(dragCurve) is str:
            values = self.readDragCurve(dragCurve)
            self.cdValues = values[0]
            self.machValues = values[1]

            for i in range(len(self.cdValues)):
                self.cdValues[i] = float(self.cdValues[i])
            for i in range(len(self.machValues)):
                self.machValues[i] = float(self.machValues[i])

        elif type(dragCurve) is list:
            self.cdValues = float(dragCurve[0])
            self.machValues = float(dragCurve[1])
        else:
            sys.exit("dragCurve variable is not given in an accepted format")

    def setFuel(self, fuelName, oxName):

        self.fuel = propellant.Propellant(fuelName)
        self.ox = propellant.Propellant(oxName)

    def performCEA(self, chamberPressure, nozzleAmbientPressure=1.01325, c='equilibrium', ct='equilibrium', sResolution=0.1, sSearch=2):

        self.chamberPressure = chamberPressure
        ceaOutput = cea.idealMixtureRatioCEA(self.fuel.ceaName, self.ox.ceaName, chamberPressure, ambientPressure=nozzleAmbientPressure, condition=c, throatCondition=c, searchResolution=sResolution, startSearch=sSearch)
        self.mixtureRatio = ceaOutput["Mixture Ratio"]
        self.specificImpulse = ceaOutput["Specific Impulse"]
        self.chamberTemperature = ceaOutput["Chamber Temperature"]
        self.expansionRatio = ceaOutput["Expansion Ratio"]
        self.exhaustGamma = ceaOutput["Gamma"]

    def calculateDeltaV(self):
        
        self.deltaV = self.specificImpulse * self.G * log(self.wetMass / self.dryMass)

    def setThrust(self, thrust):

        self.thrust = thrust

    def flightSimulation(self, timeStep=0.1, variablesToPlot="All", showPlot=False, saveFileName="Simulation.png", plotAbsoluteValues=True):
        
        self.flightSimulation = flight.Flight(self.frontalArea, self.cdValues, self.machValues)
        self.flightSimulation.flightSimulation(self.wetMass, self.propellantMass, self.thrust, self.specificImpulse, timeStep=timeStep, saveAbsoluteValues=plotAbsoluteValues)
        self.flightSimulation.plotGraph(variablesToPlot=variablesToPlot, show=showPlot, fileName=saveFileName)