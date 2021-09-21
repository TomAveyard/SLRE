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
testFuelPump = Pump(isentropicEfficiency=0.7, outletPressure=150e5)
testCoolingChannels = CoolingChannels(numberOfChannels=96, wallThickness=0.9e-3, midRibThickness=1.5e-3, channelHeight=1.25e-3, wallConductivity=365, wallRoughnessHeight=0)
testRegenerativeCooling = RegenerativeCooling(thrustChamber=testThrustChamber, coolingChannels=testCoolingChannels, solverParameters=solverParameters)
testTurbine = Turbine(isentropicEfficiency=0.7, outletPressure=testThrustChamber.injectionPressure*10**5)

testFuelLine = Line(inletState=testFuelTank.outletState, massFlowRate=testThrustChamber.fuelMassFlowRate, components=[testFuelPump, testRegenerativeCooling, testTurbine])

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
testOxPump = Pump(isentropicEfficiency=0.7, outletPressure=testThrustChamber.injectionPressure*10**5)

testOxLine = Line(inletState=testOxTank.outletState, massFlowRate=testThrustChamber.oxMassFlowRate, components=[testOxPump])

testCycle = Cycle(fuelLine=testFuelLine, oxLine=testOxLine, thrustChamber=testThrustChamber)

if plot == "cycle":
    diagramStates = [testFuelTank.outletState, testFuelPump.outletState, testRegenerativeCooling.outletState, testTurbine.outletState]
    testCycleDiagram = TSCycleDiagram(testThrustChamber.fuel.name, diagramStates)
    testCycleDiagram.drawCycle(constantPressureLines=False)