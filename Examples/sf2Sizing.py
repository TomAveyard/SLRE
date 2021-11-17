from PropTools.SubSystems.Engine.ThrustChamber.thrustChamber import ThrustChamber
import numpy as np
import matplotlib.pyplot as plt
from math import pi, sqrt

numberOfPoints = 100
contractionRatio = 4

thrusts = np.linspace(0.1, 10, numberOfPoints)
exitDiameters = np.zeros(numberOfPoints)
chamberDiameters = np.zeros(numberOfPoints)

counter = 0

for i in thrusts:

    thrustChamber = ThrustChamber(fuelName='propanol', oxName='nitrous oxide', thrust=i*1e3, chamberPressure=20, fac=True, contractionRatio=4, ambientPressure=1.01325, mixtureRatioOverride=2.5)
    exitDiameters[counter] = thrustChamber.exitDiameter*1e3

    chamberDiameter = sqrt((thrustChamber.throatArea * contractionRatio) / pi) * 2
    chamberDiameters[counter] = chamberDiameter * 1e3

    counter += 1

fig, ax = plt.subplots()
ax.plot(thrusts, exitDiameters, label="Nozzle Exit Diameter")
ax.plot(thrusts, chamberDiameters, label="Chamber Diameter (CR=4)")

ax.set_title("Variation of Diameter with Thrust")
ax.set_xlabel("Thrust [kN]")
ax.set_ylabel("Diameter [mm]")
ax.legend()

plt.show()
