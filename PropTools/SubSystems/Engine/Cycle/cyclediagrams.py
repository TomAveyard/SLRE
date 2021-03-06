import matplotlib.pyplot as plt
from PropTools.SubSystems.Engine.Propellant.propellant import Propellant
from math import ceil, floor

#TODO: Update/move to dedicated diagramming module

class TSCycleDiagram:

    def __init__(self, substance: str, states: list, SPadding: float = 100, TPadding: float = 10):

        self.substance = substance
        self.states = states

        self.TMax = 0
        self.TMin = 0
        self.SMax = 0
        self.SMin = 0

        for state in self.states:
            
            # Finds maximum temp
            if state.T > self.TMax:
                self.TMax = state.T
            # Finds minimum temp
            elif state.T < self.TMin:
                self.TMin = state.T
            # Finds maximum entropy
            if state.S > self.SMax:
                self.SMax = state.S
            # Finds minimum entropy
            elif state.S < self.SMin:
                self.TMin = state.S

        # Uses above to set temp and entropy limits of graphs
        self.TLims = [self.TMin, self.TMax + TPadding]
        self.SLims = [self.SMin, self.SMax + SPadding]

        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)

    # Draws a line of constant pressure on the diagram
    def constantPressureLine(self, pressure: float, resolution: float = 50, colour: str = 'grey', lw: float = 0.1) -> None:

        xCoords = []
        yCoords = []

        # Iterates through the entropy range of the diagram
        for i in range(round(self.SLims[0]), round(self.SLims[1]), resolution):
            xCoords.append(i)
            state = Propellant(self.substance)
            # For each entropy value at the constant pressure, we find the temp 
            state.defineState("P", pressure, "S", i)
            yCoords.append(state.T)

            lastCoord = [i, state.T]

        self.ax.plot(xCoords, yCoords, color=colour, lw=lw)
        self.ax.annotate(str(round(pressure, -4) / 1e5) + " [Bar]", xy=(lastCoord[0], lastCoord[1]), xycoords='data', xytext=(0, 5), textcoords='offset points', color='Gray', alpha=0.5, fontsize="xx-small")

    # Draws a saturation curve on the diagram
    def saturationCurve(self) -> None:

        xCoords = []
        yCoords = []
        state = Propellant(self.substance)

        # Iterates through every temp to the critical point
        for t in range(ceil(state.Tmin), floor(state.Tcrit)):
            # State is defined by the temp, and a quality of 0 to find the boundary between liquid and liquid-gas mixture states
            state.defineState("T", t, "Q", 0)
            xCoords.append(state.S)
            yCoords.append(t)

        # Plots critical point
        state.defineState("T", state.Tcrit, "P", state.Pcrit)
        xCoords.append(state.S)
        yCoords.append(state.T)

        # Iterates through every temp after the critical point
        for t in range(floor(state.Tcrit), ceil(state.Tmin), -1):
            # State is defined by the temp, and a quality of 1 to find the boundary between liquid-gas mixture and gas states
            state.defineState("T", t, "Q", 1)
            xCoords.append(state.S)
            yCoords.append(t)

        self.ax.plot(xCoords, yCoords, color='black')

    # Plots the critical point on the diagram
    def criticalPoint(self, s: float = 2) -> None:

        state = Propellant(self.substance)
        state.defineState("T", state.Tcrit, "P", state.Pcrit)

        self.ax.scatter(state.S, state.T, s=s)
        self.ax.annotate("Critical Point", xy=(state.S, state.T), xycoords='data', xytext=(0, 2), textcoords='offset points', color='Gray', alpha=0.5, fontsize="xx-small")

    # For a process that involves change in pressure, temperature, and entropy, this will plot through intermediate states in a linear fashion
    def lineInterpolation(self, tempLims: list, pressureLims: list, resolution: float = 1, colour: str = 'red', lw: float = 2) -> None:

        xCoords = []
        yCoords = []

        # Finds the "jumps" in temp we will plot over defined by the requested resolution
        tempDivisions = (round(tempLims[1]) - round(tempLims[0])) / resolution
        # Finds the "jumps" in pressure we will plot over
        pressureDivision = (pressureLims[1] - pressureLims[0]) / tempDivisions
        p = pressureLims[0]

        # Finds the states for each "jump" in pressure in temp
        for t in range(round(tempLims[0]), round(tempLims[1]), resolution):

            yCoords.append(t)

            state = Propellant(self.substance)
            state.defineState("T", t, "P", p)
            xCoords.append(state.S)

            p = p + pressureDivision

        self.ax.plot(xCoords, yCoords, color=colour, lw=lw)

    # Draws the cycle
    def plotCycle(self, constantPressureLines: bool = True, saturationCurve: bool = True, criticalPoint: bool = True):
        
        tankState = self.states[0]
        pumpState = self.states[1]
        regenState = self.states[2]
        turbineState = self.states[3]


        if saturationCurve == True:

            self.saturationCurve()

        if criticalPoint == True:

            self.criticalPoint()

        for state in self.states:

            if constantPressureLines == True:

                self.constantPressureLine(state.P)
            
            self.ax.scatter(state.S, state.T, color='black', lw=1, s=1)

        self.ax.plot([tankState.S, pumpState.S], [tankState.T, pumpState.T], color='red', lw=2)
        self.lineInterpolation([pumpState.T, regenState.T], [pumpState.P, regenState.P])
        self.ax.plot([regenState.S, turbineState.S], [regenState.T, turbineState.T], color='red', lw=2)

        self.ax.set_xlabel("Entropy [J K^-1]")
        self.ax.set_ylabel("Temperature [K]")

        self.ax.set_xlim(right=self.SLims[-1])

        plt.show()
        

            
            
            
            
            


