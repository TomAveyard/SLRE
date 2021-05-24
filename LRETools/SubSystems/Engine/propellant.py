from CoolProp.CoolProp import PropsSI

class Propellant:

    def __init__(self, name):

        self.name = name
        self.setSymbols()
        # List of properties that can be defined for a propellant
        self.properties = ["T", "P", "D", "H", "S", "Q"]
        self.Tcrit = PropsSI(self.name, "Tcrit")
        self.Tmin = PropsSI(self.name, "Tmin")
        self.Pcrit = PropsSI(self.name, "pcrit")

    # Sets all names and symbols that may be needed by various other packages
    # The propellants in this list are the propellants that are currently supported by this class
    def setSymbols(self):

        if self.name == "Methane":
            self.symbol = "CH4"
            self.ceaName = "CH4"
            self.type = "Fuel"
        
        if self.name == "Oxygen":
            self.symbol = "O2"
            self.ceaName = "LOX"
            self.type = "Ox"

    # Defines the entire state of the propellant and the point specified by 2 properties of choice
    def defineState(self, property1, property1Value, property2, property2Value):

            # Temporary dictionary to store property values
            propertiesDict = {
                "T": 0,
                "P": 0,
                "D": 0,
                "H": 0,
                "S": 0,
                "Q": 0,
            }

            # Removes stated properties from a copy of the available properties so that they can looped over for solving the rest of the properties
            propertiesToSolve = self.properties.copy()
            propertiesToSolve.remove(property1)
            propertiesToSolve.remove(property2)

            # Set the value of the stated property values
            propertiesDict[property1] = property1Value
            propertiesDict[property2] = property2Value

            # Loop over the properties that need to be solved
            for property in propertiesToSolve:

                propertiesDict[property] = PropsSI(property, property1, property1Value, property2, property2Value, self.name)

            # Save the state to the object
            self.T = propertiesDict["T"]
            self.P = propertiesDict["P"]
            self.D = propertiesDict["D"]
            self.H = propertiesDict["H"]
            self.S = propertiesDict["S"]
            self.Q = propertiesDict["Q"]
