from LRETools import vehicle
from CoolProp.CoolProp import PropsSI

fuel = "CH4"
fuelCP = "Methane"
ox = "LOX"
oxCP = "Oxygen"

karmanAlpha = vehicle.Vehicle("Karman Alpha")

karmanAlpha.setFuel(fuel, ox)
karmanAlpha.setDragCharacteristics(0.11, "dragCurve.csv")

chamberPressure = 10

karmanAlpha.performCEA(chamberPressure)

totalTankVolume = 0.15
fuelTempTank = 100
oxTempTank = 66

fuelPressureTank = 3 * 10 ** 5
oxPressureTank = 3 * 10 ** 5

fuelDensityTank = PropsSI("D", "T", fuelTempTank, "P", fuelPressureTank, fuelCP)
oxDensityTank = PropsSI("D", "T", oxTempTank, "P", oxPressureTank, oxCP)

fuelMassTank = totalTankVolume / ((karmanAlpha.mixtureRatio / oxDensityTank) + (1 / fuelDensityTank))
oxMassTank = karmanAlpha.mixtureRatio * fuelMassTank

karmanAlpha.setMasses(39, 50, fuelMassTank, oxMassTank)
karmanAlpha.setThrust(3.565 * 10 ** 3)
karmanAlpha.calculateDeltaV()

karmanAlpha.flightSimulation(timeStep=0.1, showPlot=True, variablesToPlot=["Altitude", "Velocity", "Acceleration"])
