from PropTools.SubSystems.Engine.Cycle.cycle import Cycle, Line
from PropTools.SubSystems.Engine.Cycle.pipe import Pipe
from PropTools.SubSystems.Engine.Turbopump.pump import Pump
from PropTools.SubSystems.Engine.Turbopump.turbine import Turbine
from PropTools.SubSystems.Engine.Propellant.propellant import Propellant
from PropTools.SubSystems.Engine.ThrustChamber.regenerativeCooling import RegenerativeCooling, SolverParameters
from PropTools.SubSystems.Engine.ThrustChamber.thrustChamber import ThrustChamber
from PropTools.SubSystems.Engine.ThrustChamber.regenerativeCooling import RegenerativeCooling, CoolingChannels
import matplotlib.pyplot as plt
import csv

# Define thrust chamber
sfThrustChamber = ThrustChamber(fuelName='propanol', oxName='nitrous oxide', thrust=10*10**3, chamberPressure=20, fac=True, contractionRatio=4, ambientPressure=1.01325, mixtureRatioOverride=2.75)

# Define the nozzle geometry
sfThrustChamber.getConicalNozzleGeometry(numberOfPoints=100)

# Define chamber geometry
sfThrustChamber.getChamberGeometry(lStar=1.2,
                                    contractionLength=0.075,
                                    entranceRadiusOfCurvatureFactor=0.75, 
                                    throatEntranceStartAngle=-135, 
                                    numberOfPointsConverging=100,
                                    numberOfPointsStraight=50)

# Modify parameters for the solver
solverParameters = SolverParameters(bartzEquationCoefficient=0.026, coolantSideHeatTransferCorrelation="sieder tate")

# Define the states of the propellants in the tanks
fuelTank = Propellant(sfThrustChamber.fuel.name)
fuelTank.defineState("T", 298, "P", 25*10**5)
oxTank = Propellant(sfThrustChamber.ox.name)
oxTank.defineState("T", 253, "P", 25*10**5)

# Define the components on the fuel line
fuelPipe1 = Pipe(diameter=0.01905, length=0.3, surfaceRoughness=0.0007874)
fuelCoolingChannels = CoolingChannels(numberOfChannels=90, 
                                    wallThickness=1e-3, 
                                    ribThickness=1e-3, 
                                    channelHeight=2e-3, 
                                    wallConductivity=14, 
                                    wallRoughnessHeight=6e-6,
                                    helixAngle=90)
fuelRegenCooling = RegenerativeCooling(thrustChamber=sfThrustChamber, coolingChannels=fuelCoolingChannels, solverParameters=solverParameters)
fuelPipe2 = Pipe(diameter=0.0127, length=0.1, surfaceRoughness=0.0007874)

# Define the components on the oxidiser line
oxPipe1 = Pipe(diameter=0.01905, length=0.2)

# Define the lines with above components
fuelLine = Line(inletState=fuelTank, massFlowRate=sfThrustChamber.fuelMassFlowRate, components=[fuelPipe1, fuelRegenCooling, fuelPipe2])
oxLine = Line(inletState=oxTank, massFlowRate=sfThrustChamber.oxMassFlowRate, components=[oxPipe1])

cycle = Cycle(fuelLine = fuelLine, oxLine = oxLine, thrustChamber = sfThrustChamber)

with open("heatTransferCoefficients.csv", mode='w') as file:

    fileWriter = csv.writer(file, delimiter=',')

    for i in range(len(fuelRegenCooling.gasSideHeatTransferCoefficients)):

        fileWriter.writerow([sfThrustChamber.axialCoords[i], fuelRegenCooling.gasSideHeatTransferCoefficients[i]])
