from LRETools import vehicle
from LRETools.FlightSimulation.flight import Flight
from LRETools.SubSystems.Engine.Cycle.cyclediagrams import *

expanderExample = vehicle.Vehicle("Expander Example", "Methane", "Oxygen")

# Define engine parameters
expanderExample.engine.setMass(25)
expanderExample.engine.setChamberPressure(40e5)
expanderExample.engine.setThrust(25e3)
expanderExample.engine.findIdealMixtureRatio()
expanderExample.engine.findPropellantMassFlowRates()

# Define airframe parameters
expanderExample.airframe.setMass(30)
expanderExample.airframe.setDragCharacteristics(0.11, "dragCurve.csv")
expanderExample.airframe.setTotalTankVolume(0.03)

# Define engine cycle parameters
expanderExample.engine.cycle.setFuelTank(100, 3e5)
expanderExample.engine.cycle.setOxTank(66, 3e5)
expanderExample.engine.cycle.setInjector(1e5)
expanderExample.engine.cycle.setFuelPumpEfficiency(0.7)
expanderExample.engine.cycle.setOxPumpEfficiency(0.7)
expanderExample.engine.cycle.setTurbineEfficiency(0.7)
expanderExample.engine.cycle.setRegenCoolingEstimate(147*5, 205900, expanderExample.engine.chamberTemp, 300, 4/1e3, 400, 0.2785, 9e5)

simulate = "flight"

if simulate == "cycle":
    # Draw cycle
    expanderExample.engine.cycle.solveCycle()
    cycleStates = [expanderExample.engine.cycle.fuelTank, expanderExample.engine.cycle.fuelPump, expanderExample.engine.cycle.regenCooling, expanderExample.engine.cycle.turbine]
    cycleDiagram = TSCycleDiagram("Methane", cycleStates)
    cycleDiagram.drawCycle()

elif simulate == "flight":
    # Perform flight simulation
    expanderExample.calculateMasses()
    expanderExample.calculateDeltaV()
    
    flightSim = Flight(expanderExample)
    flightSim.flightSimulation()
    flightSim.plotGraph()
