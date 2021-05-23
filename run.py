from LRETools.SubSystems.cycle import ExpanderCycle
from LRETools.SubSystems.thermodiagrams import TSCycleDiagram
from LRETools import vehicle
from CoolProp.CoolProp import PropsSI

fuel = "CH4"
fuelCP = "Methane"
ox = "LOX"
oxCP = "Oxygen"

karmanAlpha = vehicle.Vehicle("Karman Alpha")

karmanAlpha.setFuel(fuelCP, oxCP)
karmanAlpha.setDragCharacteristics(0.11, "dragCurve.csv")

chamberPressure = 40

karmanAlpha.performCEA(chamberPressure)

totalTankVolume = 0.03
fuelTempTank = 100
oxTempTank = 66

fuelPressureTank = 3 * 10 ** 5
oxPressureTank = 3 * 10 ** 5

fuelDensityTank = PropsSI("D", "T", fuelTempTank, "P", fuelPressureTank, fuelCP)
oxDensityTank = PropsSI("D", "T", oxTempTank, "P", oxPressureTank, oxCP)

fuelMassTank = totalTankVolume / ((karmanAlpha.mixtureRatio / oxDensityTank) + (1 / fuelDensityTank))
oxMassTank = karmanAlpha.mixtureRatio * fuelMassTank

karmanAlpha.setMasses(39, 50, fuelMassTank, oxMassTank)
karmanAlpha.setThrust(40 * 10 ** 3)
karmanAlpha.calculateDeltaV()

#karmanAlpha.flightSimulation(timeStep=0.1, showPlot=True, variablesToPlot=["Altitude", "Velocity", "Acceleration"])

cycle = ExpanderCycle("Methane", "Oxygen")
cycle.setFuelTank(100, 3e5)
cycle.setOxTank(66, 3e5)
cycle.setChamberConditions(karmanAlpha.chamberPressure*1e5, karmanAlpha.chamberTemperature)
cycle.setInjector(1e5)

propellantMassFlow = karmanAlpha.thrust / (karmanAlpha.specificImpulse * karmanAlpha.G)
fuelMassFlow = propellantMassFlow / karmanAlpha.mixtureRatio
oxMassFlow = fuelMassFlow * karmanAlpha.mixtureRatio

cycle.setMassFlowRate(fuelMassFlow, oxMassFlow)
cycle.setFuelPumpEfficiency(0.7)
cycle.setOxPumpEfficiency(0.7)
cycle.setTurbineEfficiency(0.7)

T_g = karmanAlpha.chamberTemperature
T_l = 300
h_g = 147 * 20
h_l = 205900
t_w = 4 / 1e3
k = 400
A = 0.2785

cycle.setRegenCoolingEstimate(T_g, T_l, h_g, t_w, k, h_l, A, 9e5)

cycle.solveCycle()

tsDiagram = TSCycleDiagram("Methane", [cycle.fuelTank, cycle.fuelPump, cycle.regenCooling, cycle.turbine])

#tsDiagram.drawCycle()
