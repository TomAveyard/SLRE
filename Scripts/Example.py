import numpy as np
from math import sqrt
from PropTools.SubSystems.Engine.ThrustChamber.regenerativeCooling import RegenerativeCooling
from PropTools.SubSystems.Engine.ThrustChamber.thrustChamber import ThrustChamber

test = ThrustChamber('ethanol', 'oxygen', 10*10**3, 25, fac=True, CR=5, ambientPressure=1)

testRC = RegenerativeCooling(test)


"""
test.getChamberGeometry(1.05, 0.05, entranceRadiusOfCurvatureFactor=0.75, throatEntranceStartAngle=-135)
test.getRaoBellNozzleGeometry(0.6)
test.plotGeometry('thrust chamber')
"""