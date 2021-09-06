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

plot = "coolant temps"

testThrustChamber = ThrustChamber('methane', 'oxygen', 7.25*10**3, 40, mixtureRatioOverride=3.16,fac=True, CR=7.5, ambientPressure=0.65)

testThrustChamber.getChamberGeometry(1.5,
                                     0.05, 
                                     entranceRadiusOfCurvatureFactor=0.75, 
                                     throatEntranceStartAngle=-135, 
                                     numberOfPointsConverging=30,
                                     numberOfPointsStraight=10)

testThrustChamber.getRaoBellNozzleGeometry(0.7, numberOfPoints=40)
testThrustChamber.getThrustChamberCoords()

print(testThrustChamber.fuelMassFlowRate)

if plot == "thrust chamber":
    fig, ax = plt.subplots()
    ax.plot(testThrustChamber.axialCoords, testThrustChamber.radialCoords)
    plt.axis("square")
    plt.show()
    exit()

fuel = Propellant(testThrustChamber.fuel.name)
ox = Propellant(testThrustChamber.ox.name)

fuel.defineState("T", 100, "P", 3*10**5)

testFuelTank = Tank(fuel)
testFuelPump = Pump(0.7, pressureRise=57*10**5)
testCoolingChannels = CoolingChannels(72, 0.001, 0.001, 0.001, 365, 6*10**(-6))
testRegenerativeCooling = RegenerativeCooling(testThrustChamber, testCoolingChannels)
testTurbine = Turbine(0.7, outletPressure=testThrustChamber.injectionPressure*10**5)

testFuelLine = Line(testFuelTank.outletState, testThrustChamber.fuelMassFlowRate, [testFuelPump, testRegenerativeCooling, testTurbine])

if plot == "heat flux":
    fig, ax = plt.subplots()
    ax.plot(testThrustChamber.axialCoords[1:-1], testRegenerativeCooling.heatFluxes[1:-1])
    plt.show()
    exit()
elif plot == "coolant temps":
    fig, ax = plt.subplots()
    ax.plot(testThrustChamber.axialCoords[1:-1], testRegenerativeCooling.coolantBulkTemps[1:-1])
    plt.show()
    exit()
elif plot == "wall temps":
    fig, ax = plt.subplots()
    ax.plot(testThrustChamber.axialCoords[1:-1], testRegenerativeCooling.gasSideWallTemps[1:-1])
    plt.show()
    exit()
elif plot == "coolant pressure":
    fig, ax = plt.subplots()
    ax.plot(testThrustChamber.axialCoords[1:-1], testRegenerativeCooling.coolantPressures[1:-1])
    plt.show()
    exit()

ox.defineState("T", 60, "P", 3*10**5)

testOxTank = Tank(ox)
testOxPump = Pump(0.7, outletPressure=testThrustChamber.injectionPressure*10**5)

testOxLine = Line(testOxTank.outletState, testThrustChamber.oxMassFlowRate, [testOxPump])

testCycle = Cycle(testFuelLine, testOxLine, testThrustChamber)

if plot == "cycle":
    diagramStates = [testFuelTank.outletState, testFuelPump.outletState, testRegenerativeCooling.outletState, testTurbine.outletState]
    testCycleDiagram = TSCycleDiagram(testThrustChamber.fuel.name, diagramStates)
    testCycleDiagram.drawCycle(constantPressureLines=False)