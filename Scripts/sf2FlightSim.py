from PropTools.SubSystems.Engine.Propellant.propellant import Propellant
from PropTools.SubSystems.Engine.ThrustChamber.regenerativeCooling import RegenerativeCooling, SolverParameters
from PropTools.SubSystems.Engine.ThrustChamber.thrustChamber import ThrustChamber
from PropTools.SubSystems.Engine.ThrustChamber.regenerativeCooling import RegenerativeCooling, CoolingChannels
from PropTools.vehicle import Vehicle
from PropTools.SubSystems.Body.airframe import Airframe
from PropTools.FlightSimulation.flight import Flight
from math import pi

sf2ThrustChamber = ThrustChamber(fuelName='ethanol', oxName='oxygen', thrust=10*10**3, chamberPressure=35, fac=True, contractionRatio=5, ambientPressure=0.7)

sf2ThrustChamber.getChamberGeometry(lStar=1.1,
                                    contractionLength=0.05, 
                                    entranceRadiusOfCurvatureFactor=0.75, 
                                    throatEntranceStartAngle=-135, 
                                    numberOfPointsConverging=100,
                                    numberOfPointsStraight=50)

sf2ThrustChamber.getRaoBellNozzleGeometry(lengthFraction=0.8, numberOfPoints=100)
sf2ThrustChamber.getThrustChamberCoords()

fuel = Propellant(sf2ThrustChamber.fuel.name)
ox = Propellant(sf2ThrustChamber.ox.name)

fuel.defineState("T", 298, "P", 3*10**5)
ox.defineState("T", 60, "P", 3*10**5)

karmanAlpha = Vehicle("Karman Alpha", "ethanol", "oxygen")

tankDiameter = 0.30
wallThickness = 0.002
outerDiameter = tankDiameter + 2 * wallThickness

karmanAlpha.airframe.setDragCharacteristics(outerDiameter, "Scripts/dragCurve.csv")

karmanAlpha.engine.thrust = sf2ThrustChamber.thrust
karmanAlpha.engine.specificImpulse = sf2ThrustChamber.specificImpulse

karmanAlpha.airframe.setMass(30)

burnTime = 46

ethanolMassRequired = sf2ThrustChamber.fuelMassFlowRate * burnTime
oxygenMassRequired = sf2ThrustChamber.oxMassFlowRate * burnTime

ethanolVolumeRequired = ethanolMassRequired / fuel.D
oxygenVolumeRequired = oxygenMassRequired / ox.D
totalVolume = ethanolVolumeRequired + oxygenVolumeRequired
tankHeight = totalVolume / (pi * ((tankDiameter / 2) ** 2))

karmanAlpha.fuelMass = ethanolMassRequired
karmanAlpha.oxMass = oxygenMassRequired
karmanAlpha.propellantMass = karmanAlpha.fuelMass + karmanAlpha.oxMass
karmanAlpha.dryMass = karmanAlpha.airframe.mass
karmanAlpha.wetMass = karmanAlpha.dryMass + karmanAlpha.propellantMass

flight = Flight(karmanAlpha)

flight.flightSimulation(timeStep=0.5)

print("---")
print("Apogee: " + str(flight.apogee))
print("Fuel Mass Required: " + str(ethanolMassRequired))
print("Oxidiser Mass Required: " + str(oxygenMassRequired))
print("Fuel Volume Required: " + str(ethanolVolumeRequired))
print("Oxidiser Volume Required: " + str(oxygenVolumeRequired))
print("Total Cylindrical Tank Height Required: " + str(tankHeight))
print("---")
