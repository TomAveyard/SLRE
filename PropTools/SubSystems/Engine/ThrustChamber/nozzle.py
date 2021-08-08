from numpy.lib import math
import PropTools.Utils.mathsUtils as mathsUtils
from math import sqrt, pi, sin, cos, tan, asin, radians, degrees
import numpy as np
import matplotlib.pyplot as plt

class Nozzle:

    def __init__(self, expansionRatio, throatRadius, entranceRadiusOfCurvatureFactor=1.5, throatEntranceInitialAngle=-135, numberOfPoints=100):

        self.expansionRatio = expansionRatio
        self.throatRadius = throatRadius
        self.entranceRadiusOfCurvatureFactor = entranceRadiusOfCurvatureFactor
        self.entranceRadiusOfCurvature = self.entranceRadiusOfCurvatureFactor * self.throatRadius
        self.throatEntranceInitialAngle = throatEntranceInitialAngle
        self.numberOfPoints = numberOfPoints

    def getConicalLength(self, divergenceHalfAngle):

        return self.throatRadius * (sqrt(self.expansionRatio) - 1) / tan(radians(divergenceHalfAngle))

    def getGeometryArrays(self, length, numberOfPoints):

        startCoord = self.entranceRadiusOfCurvature * sin(radians(self.throatEntranceInitialAngle))
        endCoord = length

        axialCoords = np.linspace(startCoord, endCoord, numberOfPoints)
        radialCoords = np.zeros(numberOfPoints)

        return axialCoords, radialCoords

    def getEntranceCoords(self, axialCoords, radialCoords):

        i = 0

        while axialCoords[i] < 0:

            radialCoords[i] = self.throatRadius + self.entranceRadiusOfCurvature - sqrt((self.entranceRadiusOfCurvature ** 2) - (axialCoords[i] ** 2))

            i += 1

        return axialCoords, radialCoords, i

class conicalNozzle(Nozzle):

    def __init__(self, expansionRatio, throatRadius, divergenceHalfAngle=15, entranceRadiusOfCurvatureFactor=1.5, throatEntranceInitialAngle=-135, numberOfPoints=100):

        super().__init__(expansionRatio, 
                        throatRadius, 
                        entranceRadiusOfCurvatureFactor=entranceRadiusOfCurvatureFactor, 
                        throatEntranceInitialAngle=throatEntranceInitialAngle, 
                        numberOfPoints=numberOfPoints
                        )

        self.divergenceHalfAngle = divergenceHalfAngle
        self.length = self.getConicalLength(self.divergenceHalfAngle)
        self.axialCoords, self.radialCoords = self.getGeometryArrays(self.length, self.numberOfPoints)
        self.getNozzleCoords()

    def getNozzleCoords(self):

        i = self.getEntranceCoords(self.axialCoords, self.radialCoords)[-1]
        
        while i < self.numberOfPoints:

            self.radialCoords[i] = self.throatRadius + self.axialCoords[i] * tan(radians(self.divergenceHalfAngle))

            i += 1


class raoBellNozzle(Nozzle):

    def __init__(self, expansionRatio, throatRadius, lengthFraction, entranceRadiusOfCurvatureFactor=1.5, throatEntranceInitialAngle=-135, numberOfPoints=100):
        super().__init__(expansionRatio, throatRadius, 
                        entranceRadiusOfCurvatureFactor=entranceRadiusOfCurvatureFactor, 
                        throatEntranceInitialAngle=throatEntranceInitialAngle, 
                        numberOfPoints=numberOfPoints
                        )

        lengthFractions = [0.6, 0.7, 0.8, 0.9, 1]

        expansionRatios = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]

        exitWallAngles = [
            [19.3, 16.8, 15.7, 14.3, 13.8, 13.2, 13.0, 12.8, 12.8, 12.5],
            [16.0, 13.8, 12.5, 11.7, 11.0, 10.5, 10.2, 10.2, 10.1, 10.0],
            [13.0, 10.7,  9.7,  8.8,  8.6,  8.5,  8.3,  8.2,  7.8,  7.7],
            [10.8,  8.2,  7.5,  7.2,  7.0,  6.8,  6.8,  6.7,  6.3,  6.2],
            [ 8.0,  6.3,  5.7,  5.3,  5.3,  5.2,  5.0,  5.0,  4.8,  4.6]
        ]

        initialWallAngles = [
            [28.7, 32.5, 34.3, 35.3, 36.2, 36.8, 37.3, 37.8, 38.2, 38.3],
            [25.7, 28.2, 33.3, 31.8, 32.8, 33.5, 34.0, 34.5, 34.7, 35.0],
            [23.0, 25.5, 27.7, 29.0, 30.0, 30.8, 31.3, 31.7, 32.0, 32.3],
            [21.5, 23.7, 26.0, 27.3, 28.7, 29.2, 29.7, 30.0, 30.3, 30.7],
            [20.0, 22.3, 24.5, 25.7, 26.5, 27.5, 28.3, 28.7, 29.2, 29.5]
        ]

        self.lengthFraction = lengthFraction
        self.length = self.getConicalLength(15) * self.lengthFraction

        self.exitWallAngle = self.getWallAngles(lengthFractions, expansionRatios, exitWallAngles)
        self.initialWallAngle = self.getWallAngles(lengthFractions, expansionRatios, initialWallAngles)

        self.throatEntranceInitialAngle = throatEntranceInitialAngle

        self.axialCoords, self.radialCoords = self.getGeometryArrays(self.length, self.numberOfPoints)
        self.getNozzleCoords()
    
    def getWallAngles(self, lengthFractions, expansionRatios, wallAngles):

        i = 0
        while lengthFractions[i] <= self.lengthFraction:
            i += 1
        tLengthFraction = mathsUtils.findRelativeDistanceBetweenNumbers(lengthFractions[i-1], lengthFractions[i], self.lengthFraction)

        j = 0
        while expansionRatios[j] <= self.expansionRatio:
            j += 1
        tExpansionRatio = mathsUtils.findRelativeDistanceBetweenNumbers(expansionRatios[i-1], expansionRatios[i], self.expansionRatio)

        if wallAngles[i][j] > wallAngles[i-1][j]:
            tLengthFraction = -tLengthFraction
        if wallAngles[i][j+1] > wallAngles[i-1][j+1]:
            tExpansionRatio = -tExpansionRatio

        wallAngle1 = wallAngles[i][j] + tLengthFraction * (abs(wallAngles[i][j] - wallAngles[i-1][j]))
        wallAngle2 = wallAngles[i][j+1] + tLengthFraction * abs(wallAngles[i][j+1] - wallAngles[i-1][j+1])
        
        if wallAngle1 > wallAngle2:
            return wallAngle2 + tExpansionRatio * abs(wallAngle2 - wallAngle1)
        else:
            return wallAngle1 + tExpansionRatio * abs(wallAngle2 - wallAngle1)

    def getNozzleCoords(self):

        i = self.getEntranceCoords(self.axialCoords, self.radialCoords)[-1]

        angle = 0

        while angle < self.initialWallAngle:

            initialRadiusOfCurvature = self.throatRadius * 0.382

            self.radialCoords[i] = self.throatRadius + initialRadiusOfCurvature - sqrt((initialRadiusOfCurvature ** 2) - (self.axialCoords[i] ** 2))

            angle = degrees(asin(self.axialCoords[i] / initialRadiusOfCurvature))

            i += 1

        bezierStart = [self.axialCoords[i-1], self.radialCoords[i-1]]
        bezierEnd = [self.length, sqrt(self.expansionRatio) * self.throatRadius]

        print(bezierEnd)

        bezierControl = mathsUtils.lineIntersection(bezierStart, tan(radians(self.initialWallAngle)), bezierEnd, tan(radians(self.exitWallAngle)))

        bezier = mathsUtils.bezierCurve([bezierStart, bezierControl, bezierEnd])

        numberOfBezierPoints = self.numberOfPoints - i

        for point in range(numberOfBezierPoints):

            t = point / numberOfBezierPoints

            coords = bezier.findPoint(t)
            self.axialCoords[i] = coords[0]
            self.radialCoords[i] = coords[1]

            i += 1


test = raoBellNozzle(10, 0.1, 0.8)

fig, ax = plt.subplots()
ax.plot(test.axialCoords, test.radialCoords)
plt.axis('square')
plt.show()
