from math import pi, sqrt, asin, tan, degrees, radians
import numpy as np
from PropTools.Utils import mathsUtils

class CombustionChamber:

    def __init__(self, 
        lStar, 
        throatRadius, 
        contractionRatio, 
        contractionLength, 
        entranceRadiusOfCurvatureFactor=1.5, 
        throatEntranceStartAngle=(-135), 
        numberOfPoints=100):

        self.lStar = lStar
        self.throatRadius = throatRadius
        self.throatArea = pi * (throatRadius ** 2)

        self.chamberVolume = self.getChamberVolume()

        self.contractionRatio = contractionRatio
        self.chamberArea = self.throatArea * contractionRatio
        self.chamberRadius = sqrt(self.chamberArea / pi)

        self.contractionLength = contractionLength

        self.entranceRadiusOfCurvatureFactor = entranceRadiusOfCurvatureFactor
        self.entranceRadiusOfCurvature = entranceRadiusOfCurvatureFactor * throatRadius

        self.throatEntranceStartAngle = throatEntranceStartAngle

        self.numberOfPoints = numberOfPoints

        self.axialCoords = np.linspace(0, -contractionLength, numberOfPoints - 1)
        self.axialCoords = np.append(self.axialCoords, 0)

        self.radialCoords = np.zeros(numberOfPoints)

        self.getChamberCoords()

        self.surfaceArea = self.getChamberSurfaceArea()


    def getChamberVolume(self):
        
        return self.lStar * self.throatArea

    # Gets the throat entrance coordinates
    def getEntranceCoords(self):

        i = 0
        angle = 0

        while angle > self.throatEntranceStartAngle:

            self.radialCoords[i] = self.throatRadius + self.entranceRadiusOfCurvature - sqrt((self.entranceRadiusOfCurvature ** 2) - (self.axialCoords[i] ** 2))

            angle = degrees(asin(self.axialCoords[i] / self.entranceRadiusOfCurvature)) - 90

            i += 1

        i -= 1
            
        return i, angle

    def getChamberCoords(self):

        i, angle = self.getEntranceCoords()

        # The following produces a quadratic bezier curve for the transition from the radial entrance section, to the cylindrical section

        bezierStart = [self.axialCoords[i], self.radialCoords[i]]
        bezierEnd = [-self.contractionLength, self.chamberRadius]

        bezierStartGradient = tan(radians(angle + 90))
        bezierEndGradient = 0

        bezierControl = mathsUtils.lineIntersection(bezierStart, bezierStartGradient, bezierEnd, bezierEndGradient)

        self.bezierPoints = [bezierStart, bezierControl, bezierEnd]

        bezier = mathsUtils.bezierCurve([bezierStart, bezierControl, bezierEnd])

        numberOfBezierPoints = self.numberOfPoints - 1 - i

        for point in range(numberOfBezierPoints):

            t = point / numberOfBezierPoints

            coords = bezier.getPoint(t)
            self.axialCoords[i] = coords[0]
            self.radialCoords[i] = coords[1]

            i += 1

        # Finds the volume of the converging section

        convergingAxialCoords = self.axialCoords[:-1]
        convergingRadialCoords = self.radialCoords[:-1]

        convergingVolume = mathsUtils.revolvedLineVolumeEstimation(convergingAxialCoords, convergingRadialCoords)

        remainingVolume = self.chamberVolume - convergingVolume

        cylindricalLength = remainingVolume / self.chamberArea

        self.axialCoords[-1] = self.axialCoords[-2] - cylindricalLength
        self.radialCoords[-1] = self.chamberRadius

        self.axialCoords = np.flip(self.axialCoords, 0)
        self.radialCoords = np.flip(self.radialCoords, 0)

    def getChamberSurfaceArea(self):

        return mathsUtils.revolvedLineSurfaceAreaEstimation(self.axialCoords, self.radialCoords)
