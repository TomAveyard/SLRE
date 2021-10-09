from copy import deepcopy
from PropTools.SubSystems.Engine.Cycle.component import Component
from PropTools.SubSystems.Engine.Propellant.propellant import Propellant
from PropTools.SubSystems.Engine.Tank.tank import Tank
from PropTools.SubSystems.Engine.Turbopump.pump import Pump
from PropTools.SubSystems.Engine.Turbopump.turbine import Turbine
from PropTools.SubSystems.Engine.ThrustChamber.thrustChamber import ThrustChamber
from PropTools.SubSystems.Engine.ThrustChamber.regenerativeCooling import RegenerativeCooling
from PropTools.SubSystems.Engine.Cycle.cyclediagrams import TSCycleDiagram

class Line:

    def __init__(self, inletState: Propellant = None, massFlowRate: float = None, components: "list[Component]" = None):

        self.inletState = inletState
        self.massFlowRate = massFlowRate
        self.components = components

        self.stateTracker = deepcopy(self.inletState)
        self.states = [deepcopy(self.stateTracker)]
        self.stationNames = ["Tank"]

        for i in components:

            if type(i) is Pump:

                i.calculate(self.stateTracker, self.massFlowRate)
                self.stationNames.append("Pump")

            elif type(i) is Turbine:

                i.calculate(self.stateTracker, self.massFlowRate)
                self.stationNames.append("Turbine")

            elif type(i) is RegenerativeCooling:

                i.calculate(self.stateTracker, self.massFlowRate)
                self.stationNames.append("Regenerative Cooling")

            self.stateTracker = deepcopy(i.outletState)
            self.states.append(deepcopy(self.stateTracker))

        self.outletState = self.states[-1]

    def outputLineStates(self, printResults=True, decimalPoints=2):

        results = """
-----------------------
{propellant} Line Results
-----------------------
Mass Flow Rate: {massFlowRate} kg/s
Inlet Temperature: {inletTemperature} K
Inlet Pressure:  {inletPressure} bar
-----------------------
        """.format(
            propellant = self.inletState.name,
            massFlowRate = round(self.massFlowRate, decimalPoints),
            inletTemperature = round(self.inletState.T, decimalPoints),
            inletPressure = round(self.inletState.P/1e5, decimalPoints)
        )

        for i in range(1, len(self.states)):

            results = results + """
{stationName} Inlet/Outlet States
-----------------------
Inlet Pressure: {inletPressure} bar
Inlet Temperature: {inletTemperature} K
Outlet Pressure: {outletPressure} bar
Outlet Temperature: {outletTemperature} K
-----------------------
            """.format(
                stationName = self.stationNames[i-1],
                inletPressure = round(self.states[i-1].P/1e5, decimalPoints),
                inletTemperature = round(self.states[i-1].T, decimalPoints),
                outletPressure = round(self.states[i].P/1e5, decimalPoints),
                outletTemperature = round(self.states[i].T, decimalPoints)
            )

        if printResults:

            print(results)

class Cycle:

    def __init__(self, fuelLine: Line = None, oxLine: Line = None, thrustChamber: ThrustChamber = None):

        self.fuelLine = fuelLine
        self.oxLine = oxLine
        self.thrustChamber = thrustChamber

        self.powerBalance = 0
        self.pumpPower = 0
        self.turbinePower = 0

        for i in fuelLine.components:

            if type(i) is Pump:

                self.pumpPower += i.power

            elif type(i) is Turbine:

                self.turbinePower += i.power

        for i in oxLine.components:

            if type(i) is Pump:

                self.pumpPower += i.power

            elif type(i) is Turbine:

                self.turbinePower += i.power

        self.powerBalance = self.turbinePower + self.pumpPower

    def plotTSDiagram(self, line: Line, SPadding: float = 100, TPadding: float = 10, includeConstantPressureLines: bool = True, includeSaturationCurve: bool = True, includeCriticalPoint: bool = True) -> None:

        plotter = TSCycleDiagram(line.inletState.name, line.states, SPadding=SPadding, TPadding=TPadding)
        plotter.plotCycle(constantPressureLines=includeConstantPressureLines, saturationCurve=includeSaturationCurve, criticalPoint=includeCriticalPoint)

        

            




        



        


