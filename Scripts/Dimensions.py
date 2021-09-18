from sys import exit
from PropTools.SubSystems.Engine.Propellant.propellant import Propellant
from PropTools.SubSystems.Engine.ThrustChamber.regenerativeCooling import RegenerativeCooling
from PropTools.SubSystems.Engine.ThrustChamber.thrustChamber import ThrustChamber
from PropTools.SubSystems.Engine.ThrustChamber.regenerativeCooling import RegenerativeCooling, CoolingChannels
from PropTools.SubSystems.Engine.Tank.tank import Tank
from PropTools.SubSystems.Engine.Turbopump.pump import Pump
from PropTools.SubSystems.Engine.Turbopump.turbine import Turbine
from PropTools.SubSystems.Engine.Cycle.cycle import Cycle, Line
from PropTools.SubSystems.Engine.Cycle.cyclediagrams import TSCycleDiagram

import matplotlib.pyplot as plt

testThrustChamber = ThrustChamber(fuelName='ethanol', oxName='oxygen', thrust=10*10**3, chamberPressure=30, fac=True, CR=5, ambientPressure=0.6)

testThrustChamber.getChamberGeometry(lStar=1.1,
                                     contractionLength=0.065, 
                                     entranceRadiusOfCurvatureFactor=0.75, 
                                     throatEntranceStartAngle=-135, 
                                     numberOfPointsConverging=50,
                                     numberOfPointsStraight=20)

testThrustChamber.getRaoBellNozzleGeometry(lengthFraction=0.8, numberOfPoints=100)
testThrustChamber.getThrustChamberCoords()

testThrustChamber.plotGeometry()
