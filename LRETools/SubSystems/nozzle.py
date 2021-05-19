from math import pi

# Class to perform quasi 1D nozzle flow calculations
class Nozzle:

    def __init__(self, gamma):
        self.gamma = gamma

        #This fraction appears in many of the equations, so for readability it has been stored in a variable
        self.commonGammaFraction = (self.gamma - 1) / 2
    
    # For the following functions, the names indicate the fraction they calculate, with underscore representing a divide

    def A_At(self, M):

        fracUpper = 1 + (self.commonGammaFraction) * (M ** 2)
        fracLower = 1 + (self.commonGammaFraction)
        exponent = (self.gamma + 1) / (2 * (self.gamma - 1))

        out = (1 / M) * ( (2 / (self.gamma + 1)) * (1 + self.commonGammaFraction * (M ** 2))) ** exponent

        return(out)

    def rho0_rho(self, M):

        exponent = 1 / (self.gamma - 1)

        out = (1 + self.commonGammaFraction * (M ** 2)) ** exponent

        return(out)

    def T0_T(self, M):

        out = 1 + self.commonGammaFraction * (M ** 2)

        return(out)

    def P0_P(self, M):

        exponent = self.gamma / (self.gamma - 1)

        out = (1 + self.commonGammaFraction * (M ** 2)) ** exponent

        return(out)

    # Inverse of the above equation, commonly used to find perfect expansion nozzle exit mach number
    def MeFromPeP0(self, Pe, P0):

        exponent = (self.gamma - 1) / self.gamma

        out = ((1 / self.commonGammaFraction) * (((P0 / Pe) ** exponent) - 1)) ** 0.5

        return(out)
    
    def P2_P1(self, M1, M2):

        upperFrac = 1 + self.gamma * (M1 ** 2)
        lowerFrac = 1 + self.gamma * (M2 ** 2)

        out = upperFrac / lowerFrac

        return(out)

    def rho2_rho1(self, M1, M2):

        upperFrac = 1 + self.commonGammaFraction * (M2 ** 2)
        lowerFrac = 1 + self.commonGammaFraction * (M1 ** 2)

        out = (M1 / M2) * ((upperFrac / lowerFrac) ** 0.5)

        return(out)

    def T2_T1(self, M1, M2):

        upperFrac = 1 + self.commonGammaFraction * (M1 ** 2)
        lowerFrac = 1 + self.commonGammaFraction * (M2 ** 2)

        out = upperFrac / lowerFrac

        return(out)

    def diameterToArea(self, d):

        out = pi * (d ** 2) / 4

        return(out)

    def areaToDiameter(self, A):

        out = ((4 * A) / pi) ** 0.5

        return(out)