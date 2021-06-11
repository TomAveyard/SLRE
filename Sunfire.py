from LRETools import vehicle
from LRETools.SubSystems.FlightSimulation.flight import Flight
import LRETools.SubSystems.Engine.Cycle.cyclediagrams as cyclediagrams

sunfire = vehicle.Vehicle("Sunfire", "Methane", "Oxygen")

# Define engine parameters
sunfire.engine.setMass(25)
sunfire.engine.setChamberPressure(40e5)
sunfire.engine.setThrust(25e3)
sunfire.engine.findIdealMixtureRatio()
sunfire.engine.findPropellantMassFlowRates()

# Define airframe parameters
sunfire.airframe.setMass(30)
sunfire.airframe.setDragCharacteristics(0.11, "dragCurve.csv")
sunfire.airframe.setTotalTankVolume(0.03)

# Define engine cycle parameters
sunfire.engine.cycle.setFuelTank(100, 3e5)
sunfire.engine.cycle.setOxTank(66, 3e5)
sunfire.engine.cycle.setInjector(1e5)
sunfire.engine.cycle.setFuelPumpEfficiency(0.7)
sunfire.engine.cycle.setOxPumpEfficiency(0.7)
sunfire.engine.cycle.setTurbineEfficiency(0.7)
sunfire.engine.cycle.setRegenCoolingEstimate(147*5, 205900, sunfire.engine.chamberTemp, 300, 4/1e3, 400, 0.2785, 9e5)

simulate = "flight"

if simulate == "cycle":
    # Draw cycle
    sunfire.engine.cycle.solveCycle()
    cycleStates = [sunfire.engine.cycle.fuelTank, sunfire.engine.cycle.fuelPump, sunfire.engine.cycle.regenCooling, sunfire.engine.cycle.turbine]
    cycleDiagram = cyclediagrams.TSCycleDiagram("Methane", cycleStates)
    cycleDiagram.drawCycle()

elif simulate == "flight":
    # Perform flight simulation
    sunfire.calculateMasses()
    sunfire.calculateDeltaV()
    
    flightSim = Flight(sunfire)
    flightSim.flightSimulation()
    flightSim.plotGraph()
