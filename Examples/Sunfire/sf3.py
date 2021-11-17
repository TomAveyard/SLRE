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
sfThrustChamber = ThrustChamber(fuelName='propanol', oxName='nitrous oxide', thrust=10*10**3, chamberPressure=20, fac=True, contractionRatio=4, ambientPressure=1.01325, mixtureRatioOverride=4)

# Define chamber geometry
sfThrustChamber.getChamberGeometry(lStar=1.2,
                                    contractionLength=0.075,
                                    entranceRadiusOfCurvatureFactor=0.75, 
                                    throatEntranceStartAngle=-135, 
                                    numberOfPointsConverging=100,
                                    numberOfPointsStraight=50)

# Define the nozzle geometry
sfThrustChamber.getConicalNozzleGeometry(numberOfPoints=100)

# Get the thrust chamber coordinates
sfThrustChamber.getThrustChamberCoords()

# Modify parameters for the solver
solverParameters = SolverParameters(bartzEquationCoefficient=0.026, coolantSideHeatTransferCorrelation="sieder tate")

# Define the states of the propellants in the tanks
fuelTank = Propellant(sfThrustChamber.fuel.name)
fuelTank.defineState("T", 298, "P", 3*10**5)
oxTank = Propellant(sfThrustChamber.ox.name)
oxTank.defineState("T", 253, "P", 25*10**5)

# Define the components on the fuel line
fuelPipe1 = Pipe(diameter=0.01905, length=0.3, surfaceRoughness=0.0007874)
fuelPump = Pump(isentropicEfficiency=0.5, outletPressure=sfThrustChamber.injectionPressure+22.5e5)
fuelPipe2 = Pipe(diameter=0.01905, length=0.25, surfaceRoughness=0.0007874)
fuelCoolingChannels = CoolingChannels(numberOfChannels=90, 
                                    wallThickness=1e-3, 
                                    ribThickness=1e-3, 
                                    channelHeight=2e-3, 
                                    wallConductivity=237, 
                                    wallRoughnessHeight=6e-6,
                                    helixAngle=90)
fuelRegenCooling = RegenerativeCooling(thrustChamber=sfThrustChamber, coolingChannels=fuelCoolingChannels, solverParameters=solverParameters)
fuelPipe3 = Pipe(diameter=0.0127, length=0.1, surfaceRoughness=0.0007874)

# Define the components on the oxidiser line
oxPipe1 = Pipe(diameter=0.01905, length=0.2)

# Define the lines with above components
fuelLine = Line(inletState=fuelTank, massFlowRate=sfThrustChamber.fuelMassFlowRate, components=[fuelPipe1, fuelPump, fuelPipe2, fuelRegenCooling, fuelPipe3])
oxLine = Line(inletState=oxTank, massFlowRate=sfThrustChamber.oxMassFlowRate, components=[oxPipe1])

cycle = Cycle(fuelLine = fuelLine, oxLine = oxLine, thrustChamber = sfThrustChamber)

# Print selected results
print("----------")
print("Engine Performance")
print("----------")
print(f"Specific Impulse: {round(sfThrustChamber.specificImpulse, 2)} s")
print("----------")
print("----------")
print("Fuel Line Component Outputs")
print("----------")
fuelPipe1.printResults(label="Fuel Pipe 1")
fuelPump.printResults(label="Fuel Pump")
fuelPipe2.printResults(label="Fuel Pipe 2")
fuelRegenCooling.printResults(label="Fuel Regen Cooling")
fuelPipe3.printResults(label="Fuel Pipe 3")
print("----------")
print("Oxidiser Line Component Outputs")
print("----------")
oxPipe1.printResults(label="Oxidiser Pipe 1")
