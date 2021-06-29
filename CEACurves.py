import PropTools.SubSystems.Engine.engine as eng
import matplotlib.pyplot as plt
import numpy as np

fuel = "Ethanol"
ox = "Hydrogen Peroxide 90"

eng = eng.Engine(fuel, ox)

chamberPressureRange = [30, 120]
numberOfPoints = 50


chamberPressuresToCheck = np.linspace(chamberPressureRange[0], chamberPressureRange[1], numberOfPoints)

idealMixtureRatios = []
chamberTemps = []
specificImpulses = []

for i in chamberPressuresToCheck:
    eng.setChamberPressure(i)
    eng.findIdealMixtureRatio(searchResolution=0.01)

    idealMixtureRatios.append(eng.idealMixtureRatio)
    chamberTemps.append(eng.chamberTemp)
    specificImpulses.append(eng.specificImpulse)

    print("Calculating @ Chamber Pressure = " + str(i))

chamberPressuresToCheck = list(chamberPressuresToCheck)

fig, axs = plt.subplots(3, 1)

axs[0].plot(chamberPressuresToCheck, idealMixtureRatios)
axs[0].set_ylabel("Ideal Mixture Ratio")

axs[1].plot(chamberPressuresToCheck, chamberTemps)
axs[1].set_ylabel("Chamber Temperature [K]")

axs[2].plot(chamberPressuresToCheck, specificImpulses)
axs[2].set_ylabel("Specific Impulse [s]")

axs[2].set_xlabel("Chamber Pressure [Bar]")

plt.suptitle("Fuel: " + fuel + " | Ox: " + ox)
plt.show()
