import sys
from math import pi, sqrt

def areaToRadius(area):

    return sqrt(area / pi)

def radiusToArea(radius):

    return pi * (radius ** 2)

def linearInterpolation(lowerPoint, upperPoint, x=None, y=None):

    if x == None and y == None:
        sys.exit("No middle point value given")
    elif x != None and y != None:
        sys.exit("Middle point is fully defined")
    
    if x != None:

        y = (((x - lowerPoint[0]) * (upperPoint[1] - lowerPoint[1])) / (upperPoint[0] - lowerPoint[0])) + lowerPoint[1]

    elif y != None:

        x = (((y - lowerPoint[1]) * (upperPoint[0] - lowerPoint[0])) / (upperPoint[1] - lowerPoint[1])) + lowerPoint[0]

    return [x, y]

def lineIntersection(point1, gradient1, point2, gradient2):

    xIntersection = (point2[1] - point1[1] - (point2[0] * gradient2) + (point1[0] * gradient1)) / (gradient1 - gradient2)
    yIntersection = (gradient1 * xIntersection) + (point1[1] - (point1[0] * gradient1))

    return [xIntersection, yIntersection]

# Generates a bezier curve given a list of control point coordinates
class bezierCurve:

    def __init__(self, controlPoints):
        
        self.controlPoints = controlPoints
        self.n = len(controlPoints) - 1

    # Binomial lookup table, more efficient than generating each time
    def binomial(self, k):

        self.pascalsTriangle = [[1],
                                [1,1],
                                [1,2,1],
                                [1,3,3,1],
                                [1,4,6,4,1],
                                [1,5,10,10,5,1],
                                [1,6,15,20,15,6,1]]

        # Generates new lines if it's not already in the above table
        while self.n > len(self.pascalsTriangle):
            nextRow = []
            nextRow.append(1)
            for i in range(len(self.pascalsTriangle[-1]) - 1):
                nextRow.append(self.pascalsTriangle[-1][i] + self.pascalsTriangle[-1][i+1])
            nextRow.append(1)
            self.pascalsTriangle.append(nextRow)

        return self.pascalsTriangle[self.n][k]

    # Finds the x/y coordinate for a given t, intended as a helper function for the findPoint function 
    def bezier(self, t, coord):

        sum = 0

        for i in range(self.n + 1):

            if coord == 'x':
                weight = self.controlPoints[i][0]
            elif coord == 'y':
                weight = self.controlPoints[i][1]
            else:
                sys.exit("Invalid coordinate type given")

            sum = sum + weight * self.binomial(i) * ((1-t) ** (self.n - i)) * (t ** i)

        return sum

    # Returns the x,y coordinates for a given t
    def getPoint(self, t):

        if t > 1 or t < 0:
            sys.exit("t has to be 0 <= t <= 1")

        x = self.bezier(t, 'x')
        y = self.bezier(t, 'y')

        return [x,y]

# Estimates the revolved volume of an abitrary line from its coordinates by adding up the volume of the truncated cone between each coordinate
def revolvedLineVolumeEstimation(axialCoords, radialCoords):

    if len(axialCoords) != len(radialCoords):

        sys.exit("Length of coordinate arrays must be equal")

    i = 0
    end = len(axialCoords) - 1
    sum = 0
    
    while i < end:

        r1 = radialCoords[i]
        r2 = radialCoords[i+1]
        l = abs(axialCoords[i+1] - axialCoords[i])

        sectionVolume = (1/3) * pi * ((r1 ** 2) + (r1 * r2) + (r2 ** 2)) * l

        sum += sectionVolume
        i += 1

    return sum

def revolvedLineSurfaceAreaEstimation(axialCoords, radialCoords):

    if len(axialCoords) != len(radialCoords):

        sys.exit("Length of coordinate arrays must be equal")

    i = 0
    end = len(axialCoords) - 1
    sum = 0

    while i < end:

        deltax = abs(axialCoords[i+1] - axialCoords[i])
        deltay = abs(radialCoords[i+1] - radialCoords[i])

        sum += sqrt((deltax ** 2) + (deltay ** 2))
        i += 1

    return sum

def distanceBetweenTwoPoints(point1, point2):

    deltax = abs(point2[0] - point1[0])
    deltay = abs(point2[1] - point1[1])

    return sqrt((deltax ** 2) + (deltay ** 2))

# Gets the distance of a number between two other numbers as a fraction of 1
# E.g. if it's 0.25 of the way between num1 and num2, then this will return 0.25
def getRelativeDistanceBetweenNumbers(num1, num2, numx):

    return (numx - num1) / (num2 - num1)

# Method of finding the radius of curvature from 3 points
# Detailed at http://paulbourke.net/geometry/circlesphere/
def radiusOfCurvature3Points2D(point1, point2, point3, returnCenterPoint=False):

    try:
        ma = (point2[1] - point1[1]) / (point2[0] - point1[0])
        mb = (point3[1] - point2[1]) / (point3[0] - point2[0])

        if (mb - ma) < 1e-12:
            return None, None

        centrex = (ma * mb * (point1[1] - point3[1]) + mb * (point1[0] + point2[0]) - ma * (point2[0] + point3[0])) / (2 * (mb - ma))
        centrey = - (1 / ma) * (centrex - ((point1[0] + point2[0]) / 2)) + ((point1[1] + point2[1]) / 2)

        radius = distanceBetweenTwoPoints([centrex, centrey], point2)

        if returnCenterPoint:
            return radius, [centrex, centrey]
        else:
            return radius

    except:
        return None, None

