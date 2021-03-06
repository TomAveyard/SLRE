from math import pi, sqrt, asin, tan, degrees, radians
import numpy as np
from PropTools.Utils import mathsUtils

class CombustionChamber:

    def __init__(self, 
        lStar: float = None, 
        throatRadius: float = None, 
        contractionRatio: float = None, 
        contractionLength: float = None, 
        entranceRadiusOfCurvatureFactor: float = 1.5, 
        throatEntranceStartAngle: float = (-135), 
        numberOfPointsConverging: int = 100,
        numberOfPointsStraight: int = 100):

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

        self.numberOfPointsConverging = numberOfPointsConverging
        self.numberOfPointsStraight = numberOfPointsStraight

        axialCoordsConverging = np.linspace(0, -contractionLength, numberOfPointsConverging)
        axialCoordsStraight = np.zeros(self.numberOfPointsStraight)
        self.axialCoords = np.append(axialCoordsConverging, axialCoordsStraight)

        self.radialCoords = np.zeros(len(self.axialCoords))

        self.getChamberCoords()

        self.surfaceArea = self.getChamberSurfaceArea()


    def getChamberVolume(self) -> None:
        
        return self.lStar * self.throatArea

    # Gets the throat entrance coordinates
    def getEntranceCoords(self) -> None:

        i = 0
        angle = 0

        while angle > self.throatEntranceStartAngle:

            self.radialCoords[i] = self.throatRadius + self.entranceRadiusOfCurvature - sqrt((self.entranceRadiusOfCurvature ** 2) - (self.axialCoords[i] ** 2))

            angle = degrees(asin(self.axialCoords[i] / self.entranceRadiusOfCurvature)) - 90

            i += 1

        i -= 1
            
        return i, angle

    def getChamberCoords(self) -> None:

        i, angle = self.getEntranceCoords()

        # The following produces a quadratic bezier curve for the transition from the radial entrance section, to the cylindrical section

        bezierStart = [self.axialCoords[i], self.radialCoords[i]]
        bezierEnd = [-self.contractionLength, self.chamberRadius]

        bezierStartGradient = tan(radians(angle + 90))
        bezierEndGradient = 0

        bezierControl = mathsUtils.lineIntersection(bezierStart, bezierStartGradient, bezierEnd, bezierEndGradient)

        self.bezierPoints = [bezierStart, bezierControl, bezierEnd]

        bezier = mathsUtils.bezierCurve([bezierStart, bezierControl, bezierEnd])

        numberOfBezierPoints = self.numberOfPointsConverging - i

        for point in range(numberOfBezierPoints):

            t = point / numberOfBezierPoints

            coords = bezier.getPoint(t)
            self.axialCoords[i] = coords[0]
            self.radialCoords[i] = coords[1]

            i += 1

        # Finds the volume of the converging section

        convergingAxialCoords = self.axialCoords[:self.numberOfPointsConverging]
        convergingRadialCoords = self.radialCoords[:self.numberOfPointsConverging]

        convergingVolume = mathsUtils.revolvedLineVolumeEstimation(convergingAxialCoords, convergingRadialCoords)

        remainingVolume = self.chamberVolume - convergingVolume

        cylindricalLength = remainingVolume / self.chamberArea

        # Add one to the numberOfPointsStraight to as there is a duplicate point where the straight and converging sections meet, which needs to be omitted
        axialCoordsStraight = np.linspace(convergingAxialCoords[-1], -cylindricalLength, self.numberOfPointsStraight + 1)

        for i in range(self.numberOfPointsStraight):

            self.axialCoords[self.numberOfPointsConverging + i] = axialCoordsStraight[i + 1]
            self.radialCoords[self.numberOfPointsConverging + i] = self.chamberRadius

        self.axialCoords = np.flip(self.axialCoords, 0)
        self.radialCoords = np.flip(self.radialCoords, 0)

    def getChamberSurfaceArea(self) -> None:

        return mathsUtils.revolvedLineSurfaceAreaEstimation(self.axialCoords, self.radialCoords)