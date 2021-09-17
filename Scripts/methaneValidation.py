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

plot = "heat flux"

testThrustChamber = ThrustChamber('methane', 'oxygen', 20.5*10**3, 56, mixtureRatioOverride=3.5, fac=True, CR=7.5, ambientPressure=0.7)

testThrustChamber.getChamberGeometry(2.25,
                                     0.065, 
                                     entranceRadiusOfCurvatureFactor=0.75, 
                                     throatEntranceStartAngle=-135, 
                                     numberOfPointsConverging=100,
                                     numberOfPointsStraight=50)

#testThrustChamber.getRaoBellNozzleGeometry(0.8, numberOfPoints=100)
testThrustChamber.getConicalNozzleGeometry(numberOfPoints=300, divergenceHalfAngle=20)
testThrustChamber.getThrustChamberCoords()

if plot == "thrust chamber":
    fig, ax = plt.subplots()
    ax.plot(testThrustChamber.axialCoords, testThrustChamber.radialCoords)
    plt.axis("square")
    plt.show()
    exit()

fuel = Propellant(testThrustChamber.fuel.name)
ox = Propellant(testThrustChamber.ox.name)

fuel.defineState("T", 108, "P", 3*10**5)

solverParameters = SolverParameters(bartzEquationCoefficient=0.023*0.33, coolantSideHeatTransferCorrelation="sieder-tate", includeRoughnessCorrection=False)

testFuelTank = Tank(fuel)
testFuelPump = Pump(0.7, outletPressure=150e5)
testCoolingChannels = CoolingChannels(96, 0.9e-3, 1.5e-3, 1.25e-3, 365, 0)
testRegenerativeCooling = RegenerativeCooling(testThrustChamber, testCoolingChannels, solverParameters)
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
elif plot == "gas side wall temps":
    fig, ax = plt.subplots()
    ax.plot(testThrustChamber.axialCoords[1:-1], testRegenerativeCooling.gasSideWallTemps[1:-1])
    plt.show()
    exit()
elif plot == "adiabatic wall temps":
    fig, ax = plt.subplots()
    ax.plot(testThrustChamber.axialCoords[1:-1], testRegenerativeCooling.adiabaticWallTemps[1:-1])
    plt.show()
    exit()
elif plot == "coolant pressure":
    fig, ax = plt.subplots()
    ax.plot(testThrustChamber.axialCoords[1:-1], testRegenerativeCooling.coolantPressures[1:-1])
    plt.show()
    exit()
elif plot == "reynold numbers":
    fig, ax = plt.subplots()
    ax.plot(testThrustChamber.axialCoords[1:-1], [i / 10**4 for i in testRegenerativeCooling.coolantReynoldsNumbers[1:-1]])
    plt.show()
    exit()
elif plot == "nusselt numbers":
    fig, ax = plt.subplots()
    ax.plot(testThrustChamber.axialCoords[1:-1], testRegenerativeCooling.coolantNusseltNumbers[1:-1])
    plt.show()
    exit()
elif plot == "prandtl numbers":
    fig, ax = plt.subplots()
    ax.plot(testThrustChamber.axialCoords[1:-1], testRegenerativeCooling.coolantPrandtlNumbers[1:-1])
    plt.show()
    exit()
elif plot == "gas side heat transfer coefficients":
    fig, ax = plt.subplots()
    ax.plot(testThrustChamber.axialCoords[1:-1], testRegenerativeCooling.gasSideHeatTransferCoefficients[1:-1])
    plt.show()
    exit()
elif plot == "coolant side heat transfer coefficients":
    fig, ax = plt.subplots()
    ax.plot(testThrustChamber.axialCoords[1:-1], testRegenerativeCooling.coolantSideHeatTransferCoefficients[1:-1])
    plt.show()
    exit()
elif plot == "channel area":
    fig, ax = plt.subplots()
    ax.plot(testThrustChamber.axialCoords[1:-1], testRegenerativeCooling.channelAreas[1:-1])
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