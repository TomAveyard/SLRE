from PropTools.SubSystems.Engine.Propellant.propellant import Propellant
from PropTools.vehicle import Vehicle
from PropTools.SubSystems.Body.airframe import Airframe
from PropTools.SubSystems.Engine.ThrustChamber.thrustChamber import ThrustChamber
from PropTools.FlightSimulation.flight import Flight
from math import pi

rocket = Vehicle("Test", "ethanol", "oxygen")

diameter = 0.3

rocket.airframe.setDragCharacteristics(0.32, "Scripts/dragCurve.csv")

thrustChamber = ThrustChamber("ethanol", 'oxygen', 10 * 10 ** 3, 30, mixtureRatioOverride=1.8, fac=True, CR=3.5)

rocket.engine.thrust = thrustChamber.thrust
rocket.engine.specificImpulse = thrustChamber.specificImpulse

rocket.airframe.setMass(60)

burnTime = 50

ethanol = Propellant("ethanol")
oxygen = Propellant("oxygen")

ethanol.defineState("T", 250, "P", 3*10**5)
oxygen.defineState("T", 100, "P", 3*10**5)

ethanolMassRequired = thrustChamber.fuelMassFlowRate * burnTime
oxygenMassRequired = thrustChamber.oxMassFlowRate * burnTime

ethanolVolumeRequired = ethanolMassRequired / ethanol.D
oxygenVolumeRequired = oxygenMassRequired / oxygen.D
totalVolume = ethanolVolumeRequired + oxygenVolumeRequired
tankHeight = totalVolume / (pi * ((diameter/2) ** 2))

rocket.fuelMass = ethanolMassRequired
rocket.oxMass = oxygenMassRequired
rocket.propellantMass = rocket.fuelMass + rocket.oxMass
rocket.dryMass = rocket.airframe.mass
rocket.wetMass = rocket.dryMass + rocket.propellantMass

flight = Flight(rocket)

flight.flightSimulation()
#flight.plotGraph(variablesToPlot=["Altitude", "Acceleration"])

#print("Apogee: " + str(flight.apogee))
#print("LiftOff TWR: " + str(rocket.engine.thrust / rocket.wetMass))
#print("Total Propellant Weight: " + str(rocket.propellantMass))
#print("Propellant Mass Flow Rate: " + str(thrustChamber.propellantMassFlowRate))
#print("Total Tank Volume: " + str(totalVolume))
#print("Total Tank Height: " + str(tankHeight))
#print("Delta V: " + str(flight.totalDragDeltaV))

print(flight.apogee)
print(ethanolMassRequired)
print(oxygenMassRequired)
