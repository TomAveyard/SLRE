from PropTools.SubSystems.Engine.Cycle.cycle import Cycle, Line
from PropTools.SubSystems.Engine.Cycle.pipe import Pipe
from PropTools.SubSystems.Engine.Turbopump.pump import Pump
from PropTools.SubSystems.Engine.Turbopump.turbine import Turbine
from PropTools.SubSystems.Engine.Propellant.propellant import Propellant
from PropTools.SubSystems.Engine.ThrustChamber.regenerativeCooling import RegenerativeCooling, SolverParameters
from PropTools.SubSystems.Engine.ThrustChamber.thrustChamber import ThrustChamber
from PropTools.SubSystems.Engine.ThrustChamber.regenerativeCooling import RegenerativeCooling, CoolingChannels
import matplotlib.pyplot as plt

# Define thrust chamber
sf2ThrustChamber = ThrustChamber(fuelName='propanol', oxName='nitrous oxide', thrust=10*10**3, chamberPressure=20, fac=True, contractionRatio=4, ambientPressure=1.01325, mixtureRatioOverride=4)

# Define chamber geometry
sf2ThrustChamber.getChamberGeometry(lStar=1.2,
                                    contractionLength=0.075,
                                    entranceRadiusOfCurvatureFactor=0.75, 
                                    throatEntranceStartAngle=-135, 
                                    numberOfPointsConverging=100,
                                    numberOfPointsStraight=50)

# Define the nozzle geometry
sf2ThrustChamber.getConicalNozzleGeometry(numberOfPoints=100)

# Get the thrust chamber coordinates
sf2ThrustChamber.getThrustChamberCoords()

# Modify parameters for the solver
solverParameters = SolverParameters(bartzEquationCoefficient=0.026, coolantSideHeatTransferCorrelation="sieder tate")

# Define the states of the propellants in the tanks
fuelTank = Propellant(sf2ThrustChamber.fuel.name)
fuelTank.defineState("T", 298, "P", 3*10**5)
oxTank = Propellant(sf2ThrustChamber.ox.name)
oxTank.defineState("T", 253, "P", 3*10**5)

# Define the components on the fuel line
fuelPipe1 = Pipe(diameter=0.01905, length=0.3, surfaceRoughness=0.0007874)
fuelPump = Pump(isentropicEfficiency=0.5, outletPressure=40e5)
fuelPipe2 = Pipe(diameter=0.01905, length=0.25, surfaceRoughness=0.0007874)
fuelCoolingChannels = CoolingChannels(numberOfChannels=90, 
                                    wallThickness=1e-3, 
                                    ribThickness=1e-3, 
                                    channelHeight=2e-3, 
                                    wallConductivity=237, 
                                    wallRoughnessHeight=6e-6,
                                    helixAngle=25)
fuelRegenCooling = RegenerativeCooling(thrustChamber=sf2ThrustChamber, coolingChannels=fuelCoolingChannels, solverParameters=solverParameters)
fuelPipe3 = Pipe(diameter=0.0127, length=0.1, surfaceRoughness=0.0007874)
fuelTurbine = Turbine(isentropicEfficiency=0.5, outletPressure=sf2ThrustChamber.injectionPressure*1e5)

# Define the components on the oxidiser line
#oxPipe1 = Pipe(diameter=0.01905, length=0.2)
oxPump = Pump(isentropicEfficiency=0.5, outletPressure=sf2ThrustChamber.injectionPressure*1e5)

# Define the lines with above components
fuelLine = Line(inletState=fuelTank, massFlowRate=sf2ThrustChamber.fuelMassFlowRate, components=[fuelPipe1, fuelPump, fuelPipe2, fuelRegenCooling, fuelPipe3, fuelTurbine])
oxLine = Line(inletState=oxTank, massFlowRate=sf2ThrustChamber.oxMassFlowRate, components=[oxPump])

cycle = Cycle(fuelLine = fuelLine, oxLine = oxLine, thrustChamber = sf2ThrustChamber)

# Print selected results
print("---")
print("Engine Performance")
print("---")
print(f"Specific Impulse: {round(sf2ThrustChamber.specificImpulse, 2)} s")
print("---")
print("Regenerative Cooling")
print("---")
print(f"Max Heat Flux: {round(max(fuelRegenCooling.heatFluxes/1e3), 2)} kW")
print(f"Max Wall Temp: {round(max(fuelRegenCooling.gasSideWallTemps), 2)} K")
print(f"Coolant Temp Rise: {round(fuelRegenCooling.coolantBulkTemps[-1] - fuelRegenCooling.coolantBulkTemps[0], 2)} K")
print(f"Coolant Outlet Temp: {round(fuelRegenCooling.coolantBulkTemps[-1], 2)} K")
print(f"Coolant Pressure Drop: {round((fuelRegenCooling.coolantPressures[-1] - fuelRegenCooling.coolantPressures[1]) / 1e5, 2)} Bar")
print("---")
print("Power Balance")
print("---")
print(f"Regenerative Cooling Heat Power: {round(fuelRegenCooling.totalHeatPower/1e3, 2)} kW")
print(f"Fuel Pump Power: {round(fuelPump.power/1e3, 2)} kW")
print(f"Oxidiser Pump Power: {round(oxPump.power/1e3, 2)} kW")
print(f"Total Pump Power: {round((fuelPump.power + oxPump.power)/1e3, 2)} kW")
print(f"Turbine Power: {round(fuelTurbine.power/1e3, 2)} kW")
print(f"Power Balance: {round((fuelTurbine.power + (fuelPump.power + oxPump.power))/1e3, 2)} kW")

# Plot
#sf2ThrustChamber.plotGeometry()
#fuelRegenCooling.plotChannels()
#fuelRegenCooling.plotHeatFlux()
#cycle.plotTSDiagram(fuelLine)