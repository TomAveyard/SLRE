from math import log10, sqrt

# Basic Equations

def massFlowRateToVelocity(massFlowRate, density, area):

        return massFlowRate / (density * area)

# Nozzle Dynamics

def getSubsonicMachNumber(localArea, throatArea, gamma, accuracy=0.001):

    # Gets local mach number upstream of the throat as it is not possible to use CEA to get this
    # Uses the Newton-Raphson method to iteratively solve for mach number from area ratio to the specified accuracy
    # Method detailed here: https://www.grc.nasa.gov/www/winddocs/utilities/b4wind_guide/mach.html

    R = (localArea / throatArea) ** 2
    E = (gamma + 1) / (gamma - 1)
    P = 2 / (gamma + 1)
    Q = 1 / E

    a = P ** (1/Q)
    r = (R - 1) / (2 * a)

    X = 1 / (1 + r + sqrt(r * (r + 2)))

    def newtonRaphson(X, P, Q):

        return (P * (X - 1)) / (1 - R * (P + Q * X) ** (-P / Q))

    Xnew = newtonRaphson(X, P, Q)

    while abs(Xnew - X) > accuracy:

        X = Xnew
        Xnew = newtonRaphson(X, P, Q)

    subsonicMachNumber = sqrt(Xnew)

    return subsonicMachNumber

def getGasTemp(area, mach, gamma, chamberTemp):

    return chamberTemp / (1 + ((gamma - 1) / 2) * (mach ** 2))

def getAdiabaticWallTemp(gasTemp, mach, gamma, prandtlNumber):

    r = prandtlNumber ** (1/3)

    return gasTemp * (1 + (r * ((gamma-1) / 2) * (mach ** 2)))

# Pipe flow

def pressureLoss(frictionFactor, length, hydraulicDiameter, density, velocity):

    return frictionFactor * (length / hydraulicDiameter) * 0.5 * density * (velocity ** 2)

def colebrookEquation(channelRoughness, hydraulicDiameter, reynoldsNumber, convergenceCriteria=0.01):

    x = -2 * log10(channelRoughness / (3.7 * hydraulicDiameter))
    xPrev = x + 10

    while abs(abs(x) - abs(xPrev)) > convergenceCriteria:

        xPrev = x
        
        x = -2 * log10((channelRoughness / (3.7 * hydraulicDiameter)) + ((2.51 * xPrev) / reynoldsNumber))

    frictionFactor = (1 / x) ** 2

    return frictionFactor

