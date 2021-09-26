from math import pi, sqrt
import sys

# Class to store information about the cooling channel design
# Channels are modelled as a sector of a annulus
class CoolingChannels:

    def __init__(self,
        numberOfChannels: int = None,
        wallThickness: float = None,
        ribThickness: float = None,
        channelHeight: float = None, 
        wallConductivity: float = None, 
        wallRoughnessHeight: float = None,
        helixAngle: float = 90):

        self.numberOfChannels = numberOfChannels
        self.wallThickness = wallThickness
        self.ribThickness = ribThickness
        self.channelHeight = channelHeight
        self.wallConductivity = wallConductivity
        self.wallRoughnessHeight = wallRoughnessHeight
        self.helixAngle = helixAngle

        self.channelInstance = ChannelDimensions(self)


# Class to store the dimensions of a single instance of the cooling channels at a certain area
class ChannelDimensions:

    def __init__(self, coolingChannels: CoolingChannels = None):

        self.coolingChannels = coolingChannels

        self.bottomWidth = None
        self.midWidth = None
        self.topWidth = None
        self.bottomRadius = None
        self.midRadius = None
        self.topRadius = None
        self.totalChannelAngle = None
        self.individualChannelAngle = None
        self.totalRibAngle = None
        self.individualRibAngle = None
        self.thrustChamberRadius = None
        self.individualChannelArea = None
        self.totalChannelArea = None
        self.sideLength = None
        self.wettedPerimeter = None
        self.hydraulicDiameter = None
        self.aspectRatio = None

    def getChannelDimensions(self, thrustChamberArea: float = None) -> None:

        self.thrustChamberRadius = sqrt(thrustChamberArea / pi)

        self.bottomRadius = self.thrustChamberRadius + self.coolingChannels.wallThickness
        self.topRadius = self.bottomRadius + self.coolingChannels.channelHeight
        self.midRadius = (self.bottomRadius + self.topRadius) / 2
        self.height = self.topRadius - self.bottomRadius

        self.totalRibAngle = (360 * self.coolingChannels.ribThickness * self.coolingChannels.numberOfChannels) / (2 * pi * self.bottomRadius)
        self.individualRibAngle = self.totalRibAngle / self.coolingChannels.numberOfChannels
        self.totalChannelAngle = 360 - self.totalRibAngle

        if self.totalChannelAngle <= 0:
            sys.exit("Error: Cooling channels are negative - aborting")

        self.individualChannelAngle = self.totalChannelAngle / self.coolingChannels.numberOfChannels

        self.bottomWidth = 2 * pi * self.bottomRadius * (self.individualChannelAngle / 360)
        self.topWidth = 2 * pi * self.topRadius * (self.individualChannelAngle / 360)
        self.midWidth = 2 * pi * self.midRadius * (self.individualChannelAngle / 360)

        self.individualChannelArea = 0.5 * self.height * (self.bottomWidth + self.topWidth)

        self.totalChannelArea = self.individualChannelArea * self.coolingChannels.numberOfChannels

        self.wettedPerimeter = self.bottomWidth + self.topWidth + (self.height * 2)
        self.hydraulicDiameter = 4 * self.individualChannelArea / self.wettedPerimeter

        self.aspectRatio = self.height / self.midWidth