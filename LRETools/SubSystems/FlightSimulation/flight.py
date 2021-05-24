import sys
from ambiance import Atmosphere
import matplotlib.pyplot as plt
from numpy import save

class Flight:

    def __init__(self, frontalArea, cdValues, machValues):

        self.G = 9.81
        self.AIR_GAMMA = 1.4
        self.AIR_R = 287

        self.frontalArea = frontalArea
        self.cdValues = cdValues
        self.machValues = machValues


    def linearInterpolation(self, xValue, x, y):

        if len(x) != len(y):
            sys.exit("Invalid: Lengths of x and y lists not equal")

        if xValue > x[-1]:
            #print("Warning: xValue greater than range of given lists, last y value has been returned")
            return(y[-1])
        
        if xValue < x[0]:
            print("Warning: xValue less than range of given lists, first y value has been returned")
            return(y[0])
        
        for i in range(len(x)):
            
            if x[i] == xValue:
                return y[i]

            elif x[i] < xValue and x[i+1] > xValue:
                x1 = x[i]
                x2 = x[i+1]
                y1 = y[i]
                y2 = y[i+1]

                yValue = ((y2 - y1) * (xValue - x1) / (x2 - x1)) + y1
                return(yValue)

    def findSoundSpeed(self, T):

        out = (self.AIR_GAMMA * self.AIR_R * T) ** 0.5

        return(out)

    def findDrag(self, mach, altitude):

        CD = self.linearInterpolation(abs(mach), self.machValues, self.cdValues)
        if altitude < 80000:
            atmosphere = Atmosphere(altitude)
            soundSpeed = self.findSoundSpeed(float(atmosphere.temperature))
            velocity = mach * soundSpeed
            drag = 0.5 * CD * float(atmosphere.density[0]) * (velocity ** 2) * self.frontalArea
        else:
            drag = 0

        return(drag)

    def flightSimulation(self, vehicleMass, propellantMass, thrust, specificImpulse, timeStep=0.1, saveAbsoluteValues=True):
        
        # Initialise temporary variables for simulation loop
        altitude = 0
        acceleration = 0
        velocity = 0
        mach = 0
        burnOut = False
        drag = 0
        time = 0
        timeStepImpulse = 0
        timeStepPropellantMassConsumption = 0

        # Initialise variables to be saved for simulation loop
        self.totalDeltaV = 0
        self.totalDragDeltaV = 0
        self.totalGravityDeltaV = 0
        self.burnTime = 0

        # Initialise lists for storing values of interest
        self.timeGraph = []
        self.altitudeGraph = []
        self.velocityGraph = []
        self.machGraph = []
        self.accelerationGraph = []
        self.vehicleMassGraph = []
        self.dragGraph = []

        print("Solving...")

        while altitude >= 0:

            print("Current Simulated Altitude: " + str(round(altitude)) + "m", end="\r")

            # Add values of previous timestep to lists
            self.timeGraph.append(time)
            self.altitudeGraph.append(altitude)
            self.vehicleMassGraph.append(vehicleMass)

            # Converts to absolute values if requested
            if saveAbsoluteValues == False:
                self.velocityGraph.append(velocity)
                self.machGraph.append(mach)
                self.accelerationGraph.append(acceleration)
                self.dragGraph.append(drag)
            elif saveAbsoluteValues == True:
                self.velocityGraph.append(abs(velocity))
                self.machGraph.append(abs(mach))
                self.accelerationGraph.append(abs(acceleration))
                self.dragGraph.append(abs(drag))

            # Calculate new altitude
            altitude = altitude + velocity * timeStep

            # Calculate drag
            drag = self.findDrag(mach, altitude)

            # Calculate force on rocket for the timestep
            if velocity >= 0:
                force = thrust - self.G * vehicleMass - drag
            elif velocity < 0:
                force = thrust - self.G * vehicleMass + drag

            # Calculate the acceleration for the timestep
            acceleration = force / vehicleMass

            # Calculate velocity for the timestep
            velocity = velocity + acceleration * timeStep

            # Calculate mach number for the timestep
            if altitude < 80000:
                atmosphere = Atmosphere(altitude)
                soundSpeed = self.findSoundSpeed(float(atmosphere.temperature))
                mach = velocity / soundSpeed
            else:
                mach = 0

            # Calculate impulse experienced during timestep, and the amount of fuel that impulse required
            timeStepImpulse = thrust * timeStep
            timeStepPropellantMassConsumption = timeStepImpulse / (specificImpulse * self.G)

            # Update masses
            propellantMass = propellantMass - timeStepPropellantMassConsumption
            vehicleMass = vehicleMass - timeStepPropellantMassConsumption

            # Update thrust according to whether fuel has run out or not
            if propellantMass <= 0 and burnOut == False:
                thrust = 0
                self.burnTime = time
                burnOut = True
            elif propellantMass > 0:
                self.totalDeltaV = self.totalDeltaV + (thrust / vehicleMass) * timeStep
                self.totalDragDeltaV = self.totalDragDeltaV + (drag / vehicleMass) * timeStep
                self.totalGravityDeltaV = self.totalGravityDeltaV + (self.G) * timeStep

            # Move onto next timestep
            time = time + timeStep

        print("\nSimulation Complete")

        self.apogee = max(self.altitudeGraph)
        self.apogeeTime = self.timeGraph[self.altitudeGraph.index(max(self.altitudeGraph))]

    def plotGraph(self, variablesToPlot="All", show=False, lw=2, fileName="Simulation.png"):

        self.fig = plt.figure()
        self.ax1 = self.fig.add_subplot(111)

        if variablesToPlot == "All":
            
            self.ax1.plot(self.timeGraph, self.altitudeGraph, linewidth=lw)
            self.ax2 = self.ax1.twinx()
            self.ax2.plot(self.timeGraph, self.velocityGraph, linewidth=lw)
            self.ax2.plot(self.timeGraph, self.accelerationGraph, linewidth=lw)
            self.ax2.plot(self.timeGraph, self.machGraph, linewidth=lw)
            self.ax2.plot(self.timeGraph, self.vehicleMassGraph, linewidth=lw)
            self.ax2.plot(self.timeGraph, self.dragGraph, linewidth=lw)
        
        elif type(variablesToPlot) is list:

            if "Altitude" in variablesToPlot:
                self.ax1.plot(self.timeGraph, self.altitudeGraph, linewidth=lw, label="Altitude", color="Black")
                self.ax2 = self.ax1.twinx()
            elif "Altitude" not in variablesToPlot:
                self.ax2 = self.ax1

            if "Velocity" in variablesToPlot:
                self.ax2.plot(self.timeGraph, self.velocityGraph, linewidth=lw, label="Velocity")
            if "Acceleration" in variablesToPlot:
                self.ax2.plot(self.timeGraph, self.accelerationGraph, linewidth=lw, label="Acceleration")
            if "Mach" in variablesToPlot:
                self.ax2.plot(self.timeGraph, self.machGraph, linewidth=lw, label="Mach")
            if "Vehicle Mass" in variablesToPlot:
                self.ax2.plot(self.timeGraph, self.vehicleMassGraph, linewidth=lw, label="Vehicle Mass")
            if "Drag" in variablesToPlot:
                self.ax2.plot(self.timeGraph, self.dragGraph, linewidth=lw, label="Drag")
            
        plt.axvline(x=self.burnTime, lw=1, alpha=0.5, color="Gray")
        self.ax1.annotate("Burn Out", xy=(self.burnTime, 0), xycoords='data', xytext=(5, 5), textcoords='offset points', color='Gray', alpha=0.5)
        plt.axvline(x=self.apogeeTime, lw=1, alpha=0.5, color="Gray")
        self.ax1.annotate("Apogee", xy=(self.apogeeTime, 0), xycoords='data', xytext=(5, 5), textcoords='offset points', color='Gray', alpha=0.5)

        self.ax1.set_ylim(bottom=0)
        self.ax1.set_xlim([0, self.timeGraph[-1]])
        self.ax2.set_ylim(bottom=0)

        plt.savefig(fileName)
        
        if show == True:
            plt.show()
