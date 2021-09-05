from PropTools.SubSystems.Engine.Propellant.propellant import Propellant
from PropTools.SubSystems.Engine.ThrustChamber.regenerativeCooling import RegenerativeCooling
from PropTools.SubSystems.Engine.ThrustChamber.thrustChamber import ThrustChamber
from PropTools.SubSystems.Engine.ThrustChamber.regenerativeCooling import RegenerativeCooling, CoolingChannels

import matplotlib.pyplot as plt

testThrustChamber = ThrustChamber('ethanol', 'oxygen', 10*10**3, 30, fac=True, CR=4, ambientPressure=0.65)

testThrustChamber.getChamberGeometry(1.05,
                                     0.1, 
                                     entranceRadiusOfCurvatureFactor=0.75, 
                                     throatEntranceStartAngle=-135, 
                                     numberOfPointsConverging=200,
                                     numberOfPointsStraight=10)

testThrustChamber.getRaoBellNozzleGeometry(0.6, numberOfPoints=200)
testThrustChamber.getThrustChamberCoords()

testCoolingChannels = CoolingChannels(40, 0.005, 0.005, 0.02, 300, 5*10**(-6))

inlet = Propellant(testThrustChamber.fuel.name)
inlet.defineState("T", 298, "P", 50*10**5)

testRegenCooling = RegenerativeCooling(testThrustChamber, testCoolingChannels)
testRegenCooling.calculate(inlet, testThrustChamber.fuelMassFlowRate, convergenceCriteria=0.01)

fig, ax = plt.subplots()
ax.plot(testThrustChamber.axialCoords[1:-1], testRegenCooling.heatFluxes[1:-1])
plt.show()
