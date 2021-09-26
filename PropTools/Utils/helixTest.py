import matplotlib.pyplot as plt
from matplotlib import cm
from math import degrees, sin, cos, tan, sqrt, radians, pi
import numpy as np
from PropTools.SubSystems.Engine.ThrustChamber.thrustChamber import ThrustChamber
from PropTools.Utils.mathsUtils import rotatePoint

class Helix:

    def __init__(self, axialCoords:list = None, radialCoords: list = None, helixAngle: float = 15, startingAngle: float = 0):
        
        self.axialCoords = axialCoords
        self.radialCoords = radialCoords
        self.helixAngle = helixAngle
        self.startingAngle = startingAngle

        theta = radians(self.startingAngle)
        self.helixAngleRadians = radians(self.helixAngle)

        self.numberOfPoints = len(self.axialCoords)

        self.xCoords = np.zeros(self.numberOfPoints)
        self.yCoords = np.zeros(self.numberOfPoints)
        self.zCoords = np.zeros(self.numberOfPoints)

        # Get first point so we can access a zCoordPrev in the loop
        self.xCoords[0] = radialCoords[0] * cos(theta)
        self.yCoords[0] = radialCoords[0] * sin(theta)
        self.zCoords[0] = axialCoords[0]

        i = 1

        while i <= self.numberOfPoints - 1:

            zCoord = axialCoords[i]
            zCoordPrev = axialCoords[i-1]

            theta += ((zCoord - zCoordPrev) * pi) / ((radialCoords[i] * tan(self.helixAngleRadians)))

            xCoord = radialCoords[i] * cos(theta)
            yCoord = radialCoords[i] * sin(theta)

            self.xCoords[i] = xCoord
            self.yCoords[i] = yCoord
            self.zCoords[i] = zCoord

            i += 1
        
    def plotHelix(self, include2DContour: bool = False, numberOfChannels: int = 1, lw: float = 1):

        ax = plt.figure().add_subplot(projection='3d')

        ax.plot(self.xCoords, self.yCoords, self.zCoords, lw=lw)

        if include2DContour:
            ax.plot(self.radialCoords, self.axialCoords, zs=0, zdir='y')

        if numberOfChannels > 1:

            thetaIncrement = 360 / numberOfChannels
            theta = self.startingAngle + thetaIncrement

            for i in range(1, numberOfChannels):
  
                newxCoords = np.zeros(self.numberOfPoints)
                newyCoords = np.zeros(self.numberOfPoints)

                for j in range(self.numberOfPoints):

                    xCoord, yCoord = rotatePoint([self.xCoords[j], self.yCoords[j]], [0,0], angle=theta)
                    newxCoords[j] = xCoord
                    newyCoords[j] = yCoord
                
                ax.plot(newxCoords, newyCoords, self.zCoords, lw=lw)

                theta += thetaIncrement

        axialLims = [min(self.axialCoords), max(self.axialCoords)]
        ax.axes.set_xlim3d(left=axialLims[0], right=axialLims[1]) 
        ax.axes.set_ylim3d(bottom=axialLims[0], top=axialLims[1]) 
        ax.axes.set_zlim3d(bottom=axialLims[0], top=axialLims[1]) 
        plt.show()

testChamber = ThrustChamber(fuelName="methane", oxName="oxygen", thrust=10*10**3, chamberPressure=35, mixtureRatioOverride=3.5, fac=True, contractionRatio=5, ambientPressure=0.2)
testChamber.getRaoBellNozzleGeometry(lengthFraction=0.8, numberOfPoints=1000)
testChamber.getChamberGeometry(lStar=1.1, contractionLength=0.1, numberOfPointsConverging=300, numberOfPointsStraight=500)
testChamber.getThrustChamberCoords()

test = Helix(axialCoords=testChamber.axialCoords, radialCoords=testChamber.radialCoords, helixAngle=30)

test.plotHelix(numberOfChannels=96)






        

