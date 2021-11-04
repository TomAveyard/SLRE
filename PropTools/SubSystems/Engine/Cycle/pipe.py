from PropTools.SubSystems.Engine.Cycle.component import Component
from PropTools.SubSystems.Engine.Propellant.propellant import Propellant
from PropTools.Thermo.fluidDynamics import colebrookEquation, pressureLoss
from math import pi
import sys

class Pipe(Component):

    def __init__(self, diameter: float = None, length: float = None, surfaceRoughness: float = 6e-6, bends: list = None, convergenceCriteria: float = 0.01):

        super().__init__()
        self.diameter = diameter
        self.length = length
        self.surfaceRoughness = surfaceRoughness
        self.bends = bends
        self.convergenceCriteria = convergenceCriteria

        self.area = pi * (self.diameter / 2) ** 2

        self.inletState = None
        self.massFlowRate = None
        self.velocity = None
        self.reynoldsNumber = None
        self.frictionFactor = None
        self.pressureLoss = None
        self.outletPressure = None
        self.outletState = None

    def calculate(self, inletState: Propellant = None, massFlowRate: float = None) -> None:
        
        self.inletState = inletState
        self.massFlowRate = massFlowRate

        self.velocity = self.massFlowRate / (self.inletState.D * self.area)
        self.reynoldsNumber = (self.inletState.D * self.velocity * self.diameter) / self.inletState.viscosity

        self.frictionFactor = colebrookEquation(self.surfaceRoughness, self.diameter, self.reynoldsNumber, convergenceCriteria=self.convergenceCriteria)
        self.pressureLoss = pressureLoss(self.frictionFactor, self.length, self.diameter, self.inletState.D, self.velocity)

        # Calculate pressure losses due to bends if specified
        if self.bends != None:

            pressureLossCoefficients = []
            
            for i in self.bends:

                if i.lower() == "90 sharp" or i.lower() == "90s" or i.lower() == "90 s":

                    pressureLossCoefficients.append(1)
                
                elif i.lower() == "45 sharp" or i.lower() == "45s" or i.lower() == "45 s":

                    pressureLossCoefficients.append(0.3)

                elif i.lower() == "90 gradual" or i.lower() == "90g" or i.lower() == "90 g":

                    pressureLossCoefficients.append(0.2)

                elif i.lower() == "45 gradual" or i.lower() == "45g" or i.lower() == "45 g":

                    pressureLossCoefficients.append(0.1)

                else:

                    sys.exit("Please enter a valid bend type")

            for k in pressureLossCoefficients:

                self.pressureLoss += k * (self.inletState.D) / 2 * (self.velocity ** 2)

        self.outletPressure = self.inletState.P - self.pressureLoss

        self.outletState = Propellant(self.inletState.name)
        self.outletState.defineState("T", self.inletState.T, "P", self.outletPressure)
