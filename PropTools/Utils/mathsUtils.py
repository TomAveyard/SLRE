import sys

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

class bezierCurve:

    def __init__(self, controlPoints):
        
        self.controlPoints = controlPoints
        self.n = len(controlPoints) - 1

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


    def findPoint(self, t):

        if t > 1 or t < 0:
            sys.exit("t has to be 0 <= t <= 1")

        x = self.bezier(t, 'x')
        y = self.bezier(t, 'y')

        return [x,y]

def findRelativeDistanceBetweenNumbers(num1, num2, numx):

    return (numx - num1) / (num2 - num1)
