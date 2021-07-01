import sys
from CoolProp.CoolProp import PropsSI
import PropTools.SubSystems.Engine.Cycle.Propellant.thermoWrapper as tw
import json
from PropTools.Utils.fileHandling import openFromRelativePath

class Propellant:

    def __init__(self, name, library=None):

        self.name = name
        self.symbol = None
        self.ceaName = None
        self.type = None
        self.mixture = None
        self.components = None
        self.setSymbols()

        if library != None:
            self.library = library

        # List of properties that can be defined for a propellant
        self.properties = ["T", "P", "D", "H", "S", "Q"]

        if self.library == "CoolProp":

            self.Tcrit = PropsSI(self.name, "Tcrit")
            self.Tmin = PropsSI(self.name, "Tmin")
            self.Pcrit = PropsSI(self.name, "pcrit")

        elif self.library == "Thermo":

            if self.mixture == "False":

                self.chemical = tw.thermoWrapper(self.name)

            elif self.mixture == "True":

                substances = list(self.components.keys())
                massFractions = []

                for s in substances:
                    massFractions.append(self.components[s])
                     
                self.chemical = tw.thermoWrapper(substances, massFractions)
                
            self.chemical.calculate(T=298, P=101325)
    
            self.Tcrit = self.chemical.constants.Tcs[0]
            self.Tmin = 20 # Thermo allows temps down to 0.01K, however temps below 20K will never be realistically be required for this project
            self.Pcrit = self.chemical.constants.Pcs[0]

        else:
            sys.exit("Specified library is either not valid, or not supported")

        # Initialise variables that will be defined using other functions
        self.T = None
        self.P = None
        self.D = None
        self.H = None
        self.S = None
        self.Q = None

    # Loads settings for propellant from propellantData.json file where these are stored
    def setSymbols(self):
        
        f = openFromRelativePath("PropTools/SubSystems/Engine/Cycle/Propellant/propellantData.json")
        data = json.load(f)

        data = data[self.name.lower()]
        
        self.symbol = data["symbol"]
        self.ceaName = data["ceaName"]
        self.type = data["type"]
        self.mixture = data["mixture"]
        self.components = data["components"]
        self.library = data["defaultLibrary"]

        f.close()

    # Defines the entire state of the propellant and the point specified by 2 properties of choice
    def defineState(self, property1, property1Value, property2, property2Value):

            # Temporary dictionary to store property values
            propertiesDict = {
                "T": None,
                "P": None,
                "D": None,
                "H": None,
                "S": None,
                "Q": None,
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
                
                if self.library == "CoolProp":
                    propertiesDict[property] = PropsSI(property, property1, property1Value, property2, property2Value, self.name)

                elif self.library == "Thermo":
                    
                    self.chemical.calculate(T=propertiesDict["T"], P=propertiesDict["P"], D=propertiesDict["D"], H=propertiesDict["H"], S=propertiesDict["S"], Q=propertiesDict["Q"])
                    propertiesDict["T"] = self.chemical.T
                    propertiesDict["P"] = self.chemical.P
                    propertiesDict["D"] = self.chemical.D
                    propertiesDict["H"] = self.chemical.H
                    propertiesDict["S"] = self.chemical.S
                    propertiesDict["Q"] = self.chemical.Q
                    
            # Save the state to the object
            self.T = propertiesDict["T"]
            self.P = propertiesDict["P"]
            self.D = propertiesDict["D"]
            self.H = propertiesDict["H"]
            self.S = propertiesDict["S"]
            self.Q = propertiesDict["Q"]