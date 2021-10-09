from sys import exit
from PropTools.SubSystems.Engine.Propellant.propellant import Propellant
from PropTools.SubSystems.Engine.ThrustChamber.regenerativeCooling import RegenerativeCooling, SolverParameters
from PropTools.SubSystems.Engine.ThrustChamber.thrustChamber import ThrustChamber
from PropTools.SubSystems.Engine.ThrustChamber.regenerativeCooling import RegenerativeCooling, CoolingChannels
from PropTools.SubSystems.Engine.Tank.tank import Tank
from PropTools.SubSystems.Engine.Turbopump.pump import Pump
from PropTools.SubSystems.Engine.Turbopump.turbine import Turbine 
from PropTools.SubSystems.Engine.Cycle.cycle import Cycle, Line
from PropTools.SubSystems.Engine.Cycle.cyclediagrams import TSCycleDiagram

import matplotlib.pyplot as plt

testThrustChamber = ThrustChamber(fuelName='methane', oxName='oxygen', thrust=20.5*10**3, chamberPressure=56, mixtureRatioOverride=3.5, fac=True, contractionRatio=7.5, ambientPressure=0.7)

testThrustChamber.getChamberGeometry(lStar=2.25,
                                     contractionLength=0.065, 
                                     entranceRadiusOfCurvatureFactor=0.75, 
                                     throatEntranceStartAngle=-135, 
                                     numberOfPointsConverging=100,
                                     numberOfPointsStraight=50)

#testThrustChamber.getRaoBellNozzleGeometry(0.8, numberOfPoints=100)
testThrustChamber.getConicalNozzleGeometry(numberOfPoints=300, divergenceHalfAngle=20)
testThrustChamber.getThrustChamberCoords()

#testThrustChamber.plotGeometry()

fuel = Propellant(testThrustChamber.fuel.name)
ox = Propellant(testThrustChamber.ox.name)

fuel.defineState("T", 108, "P", 3*10**5)

solverParameters = SolverParameters(bartzEquationCoefficient=0.026, coolantSideHeatTransferCorrelation="dittus-boelter", includeRoughnessCorrection=False)

testFuelTank = Tank(fuel)
testFuelPump = Pump(isentropicEfficiency=0.7, outletPressure=150e5)
testCoolingChannels = CoolingChannels(numberOfChannels=96, wallThickness=0.9e-3, ribThickness=1.5e-3, channelHeight=1.25e-3, wallConductivity=365, wallRoughnessHeight=0, helixAngle=90)
testRegenerativeCooling = RegenerativeCooling(thrustChamber=testThrustChamber, coolingChannels=testCoolingChannels, solverParameters=solverParameters)
testTurbine = Turbine(isentropicEfficiency=0.7, outletPressure=testThrustChamber.injectionPressure*10**5)

testFuelLine = Line(inletState=testFuelTank.outletState, massFlowRate=testThrustChamber.fuelMassFlowRate, components=[testFuelPump, testRegenerativeCooling, testTurbine])

testRegenerativeCooling.plotHeatFlux()

ox.defineState("T", 60, "P", 3*10**5)

testOxTank = Tank(ox)
testOxPump = Pump(isentropicEfficiency=0.7, outletPressure=testThrustChamber.injectionPressure*10**5)

testOxLine = Line(inletState=testOxTank.outletState, massFlowRate=testThrustChamber.oxMassFlowRate, components=[testOxPump])

testCycle = Cycle(fuelLine=testFuelLine, oxLine=testOxLine, thrustChamber=testThrustChamber)

"""
diagramStates = [testFuelTank.outletState, testFuelPump.outletState, testRegenerativeCooling.outletState, testTurbine.outletState]
testCycleDiagram = TSCycleDiagram(testThrustChamber.fuel.name, diagramStates)
testCycleDiagram.drawCycle(constantPressureLines=False)
"""