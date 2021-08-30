from PropTools.SubSystems.Engine.Propellant.propellant import Propellant
from PropTools.SubSystems.Engine.ThrustChamber.regenerativeCooling import RegenerativeCooling
from PropTools.SubSystems.Engine.ThrustChamber.thrustChamber import ThrustChamber
from PropTools.SubSystems.Engine.ThrustChamber.regenerativeCooling import RegenerativeCooling, CoolingChannels

import matplotlib.pyplot as plt

testThrustChamber = ThrustChamber('ethanol', 'oxygen', 10*10**3, 35, fac=True, CR=5, ambientPressure=1)

testThrustChamber.getChamberGeometry(1.05,
                                     0.05, 
                                     entranceRadiusOfCurvatureFactor=0.75, 
                                     throatEntranceStartAngle=-135, 
                                     numberOfPointsConverging=200,
                                     numberOfPointsStraight=20)

testThrustChamber.getRaoBellNozzleGeometry(0.6, numberOfPoints=200)
testThrustChamber.getThrustChamberCoords()

testCoolingChannels = CoolingChannels(testThrustChamber.fuelMassFlowRate, 40, 0.005, 0.005, 0.02, 300, 5*10**(-6))

inlet = Propellant(testThrustChamber.fuel.name)
inlet.defineState("T", 298, "P", 50*10**5)

testRegenCooling = RegenerativeCooling(testThrustChamber, testCoolingChannels, inlet)
testRegenCooling.calculate(convergenceCriteria=0.01)

print("Pressure Loss: " + str(testRegenCooling))

fig, ax = plt.subplots()
ax.plot(testThrustChamber.axialCoords[1:-1], testRegenCooling.gasSideWallTemps[1:-1])
plt.show()