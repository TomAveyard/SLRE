from LRETools import vehicle

karmanAlpha = vehicle.Vehicle("Karman Alpha", "Methane", "Oxygen")

karmanAlpha.engine.setMass(25)
karmanAlpha.engine.setChamberPressure(40e5)
karmanAlpha.engine.setThrust(40e3)
karmanAlpha.engine.findIdealMixtureRatio()
karmanAlpha.engine.findPropellantMassFlowRates()

karmanAlpha.airframe.setMass(30)
karmanAlpha.airframe.setDragCharacteristics(0.11, "dragCurve.csv")
karmanAlpha.airframe.setTotalTankVolume(0.03)

karmanAlpha.engine.cycle.setFuelTank(100, 3e5)
karmanAlpha.engine.cycle.setOxTank(66, 3e5)

karmanAlpha.engine.cycle.setInjector(1e5)
karmanAlpha.engine.cycle.setFuelPumpEfficiency(0.7)
karmanAlpha.engine.cycle.setOxPumpEfficiency(0.7)
karmanAlpha.engine.cycle.setTurbineEfficiency(0.7)
karmanAlpha.engine.cycle.setRegenCoolingEstimate(karmanAlpha.engine.chamberTemp, 300, 147 * 20, 4/1e3, 400, 205900, 0.2785, 9e5)
karmanAlpha.engine.cycle.solveCycle()

karmanAlpha.calculateMasses()
karmanAlpha.calculateDeltaV()
karmanAlpha.flightSimulation(showPlot=True)




