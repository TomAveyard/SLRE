# Defines common parametes and methods for line components
class Component:

    def __init__(self):

        self.inletState = None
        self.outletState = None
        self.type = None

        self.printSeperator = "----------"

    def calculate(self):

        pass

    def printResults(self):

        pass