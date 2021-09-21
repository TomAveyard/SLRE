from PropTools.SubSystems.Engine.Cycle.cycle import Line
from PropTools.SubSystems.Engine.Turbopump.pump import Pump
from PropTools.SubSystems.Engine.Turbopump.turbine import Turbine
from PropTools.SubSystems.Engine.Propellant.propellant import Propellant
from PropTools.SubSystems.Engine.ThrustChamber.regenerativeCooling import RegenerativeCooling, SolverParameters
from PropTools.SubSystems.Engine.ThrustChamber.thrustChamber import ThrustChamber
from PropTools.SubSystems.Engine.ThrustChamber.regenerativeCooling import RegenerativeCooling, CoolingChannels
import matplotlib.pyplot as plt

# Modify parameters for the solver
solverParameters = SolverParameters(bartzEquationCoefficient=0.026*0.25, coolantSideHeatTransferCorrelation="sieder tate")

# Define thrust chamber
thrustChamber = ThrustChamber(fuelName='hydrogen', oxName='oxygen', thrust=1007.7*10**3, chamberPressure=100, mixtureRatioOverride=5.6, fac=True, contractionRatio=2.1, ambientPressure=4)

# Define chamber geometry
thrustChamber.getChamberGeometry(lStar=1.35,
                                contractionLength=0.35, 
                                entranceRadiusOfCurvatureFactor=0.1, 
                                throatEntranceStartAngle=-135, 
                                numberOfPointsConverging=100,
                                numberOfPointsStraight=50)

# Define the nozzle geometry
thrustChamber.getRaoBellNozzleGeometry(lengthFraction=0.6)

# Get the thrust chamber coordinates
thrustChamber.getThrustChamberCoords()

# Define the states of the propellants in the tanks
fuelTank = Propellant(thrustChamber.fuel.name)
fuelTank.defineState("T", 16, "P", 3*10**5)
oxTank = Propellant(thrustChamber.ox.name)
oxTank.defineState("T", 60, "P", 3*10**5)

# Define the components on the fuel line
fuelPump = Pump(isentropicEfficiency=0.7, outletPressure=165e5)
fuelCoolingChannels = CoolingChannels(numberOfChannels=288, wallThickness=1.7e-3, midRibThickness=1.7e-3*0.5, channelHeight=9.1e-3, wallConductivity=390, wallRoughnessHeight=1e-6)
fuelRegenCooling = RegenerativeCooling(thrustChamber=thrustChamber, coolingChannels=fuelCoolingChannels, solverParameters=solverParameters)
fuelTurbine = Turbine(isentropicEfficiency=0.7, outletPressure=thrustChamber.injectionPressure*1e5)

# Define the components on the oxidiser line
oxPump = Pump(isentropicEfficiency=0.7, outletPressure=thrustChamber.injectionPressure*1e5)

# Define the lines with above components
fuelLine = Line(inletState=fuelTank, massFlowRate=thrustChamber.fuelMassFlowRate, components=[fuelPump, fuelRegenCooling, fuelTurbine])
oxLine = Line(inletState=oxTank, massFlowRate=thrustChamber.oxMassFlowRate, components=[oxPump])

"""
# Print selected results
print("---")
print("Engine Performance")
print("---")
print(f"Specific Impulse: {round(thrustChamber.specificImpulse, 2)} s")
print("---")
print("Regenerative Cooling")
print("---")
print(f"Max Heat Flux: {round(max(fuelRegenCooling.heatFluxes/1e3), 2)} kW")
print(f"Max Wall Temp: {round(max(fuelRegenCooling.gasSideWallTemps), 2)} K")
print(f"Coolant Temp Rise: {round(fuelRegenCooling.coolantBulkTemps[-1] - fuelRegenCooling.coolantBulkTemps[1], 2)} K")
print(f"Coolant Outlet Temp: {round(fuelRegenCooling.coolantBulkTemps[-1], 2)} K")
print(f"Coolant Pressure Drop: {round((fuelRegenCooling.coolantPressures[-1] - fuelRegenCooling.coolantPressures[0]) / 1e5, 2)} Bar")
print("---")
print("Power Balance")
print("---")
print(f"Regenerative Cooling Heat Power: {round(fuelRegenCooling.totalHeatPower/1e3, 2)} kW")
print(f"Fuel Pump Power: {round(fuelPump.power/1e3, 2)} kW")
print(f"Oxidiser Pump Power: {round(oxPump.power/1e3, 2)} kW")
print(f"Total Pump Power: {round((fuelPump.power + oxPump.power)/1e3, 2)} kW")
print(f"Turbine Power: {round(fuelTurbine.power/1e3, 2)} kW")
print(f"Power Balance: {round((fuelTurbine.power - (fuelPump.power + oxPump.power))/1e3, 2)} kW")
"""

# Plot
fig, ax = plt.subplots()
ax.plot(thrustChamber.axialCoords[1:-1], fuelRegenCooling.heatFluxes[1:-1])
plt.show()