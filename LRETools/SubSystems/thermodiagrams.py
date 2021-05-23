import matplotlib.pyplot as plt
from LRETools.SubSystems.propellant import Propellant
from math import ceil, floor

class TSCycleDiagram:

    def __init__(self, substance, states, SPadding=100, TPadding=10):

        self.substance = substance
        self.states = states

        self.TMax = 0
        self.TMin = 0
        self.SMax = 0
        self.SMin = 0

        for state in self.states:
            
            if state.T > self.TMax:
                self.TMax = state.T
            elif state.T < self.TMin:
                self.TMin = state.T
            if state.S > self.SMax:
                self.SMax = state.S
            elif state.S < self.SMin:
                self.TMin = state.S

        self.TLims = [self.TMin, self.TMax + TPadding]
        self.SLims = [self.SMin, self.SMax + SPadding]

        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)

    def constantPressureLine(self, pressure, resolution=50, colour='grey', lw=0.1):

        xCoords = []
        yCoords = []

        for i in range(round(self.SLims[0]), round(self.SLims[1]), resolution):
            xCoords.append(i)
            state = Propellant(self.substance)
            state.defineState("P", pressure, "S", i)
            yCoords.append(state.T)

        self.ax.plot(xCoords, yCoords, color=colour, lw=lw)

    def saturationCurve(self):

        xCoords = []
        yCoords = []
        state = Propellant(self.substance)

        for t in range(ceil(state.Tmin), floor(state.Tcrit)):
            state.defineState("T", t, "Q", 0)
            xCoords.append(state.S)
            yCoords.append(t)

        state.defineState("T", state.Tcrit, "P", state.Pcrit)
        xCoords.append(state.S)
        yCoords.append(state.T)

        for t in range(floor(state.Tcrit), ceil(state.Tmin), -1):
            state.defineState("T", t, "Q", 1)
            xCoords.append(state.S)
            yCoords.append(t)

        self.ax.plot(xCoords, yCoords, color='black')

    def criticalPoint(self):

        state = Propellant(self.substance)
        state.defineState("T", state.Tcrit, "P", state.Pcrit)

        self.ax.scatter(state.S, state.T)

    def lineInterpolation(self, tempLims, pressureLims, resolution=1, colour='purple', lw=2):

        xCoords = []
        yCoords = []

        tempDivisions = (round(tempLims[1]) - round(tempLims[0])) / resolution
        pressureDivision = (pressureLims[1] - pressureLims[0]) / tempDivisions
        p = pressureLims[0]

        for t in range(round(tempLims[0]), round(tempLims[1]), resolution):

            yCoords.append(t)

            state = Propellant(self.substance)
            state.defineState("T", t, "P", p)
            xCoords.append(state.S)

            p = p + pressureDivision

        self.ax.plot(xCoords, yCoords, color=colour, lw=lw)

    def drawCycle(self, constantPressureLines=True, saturationCurve=True, criticalPoint=True):
        
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
            
            self.ax.scatter(state.S, state.T, color='black', lw=1)

        self.ax.plot([tankState.S, pumpState.S], [tankState.T, pumpState.T], color='red', lw=2)
        self.lineInterpolation([pumpState.T, regenState.T], [pumpState.P, regenState.P])
        self.ax.plot([regenState.S, turbineState.S], [regenState.T, turbineState.T], color='red', lw=2)

        plt.show()
        

            
            
            
            
            


