from posixpath import relpath
import sys
import csv
from PropTools.Utils.fileHandling import openFromRelativePath

class Airframe:

    def __init__(self):
        
        # Initialise variables that will be defined using other functions
        self.mass = None
        self.frontalArea = None
        self.cdValues = None
        self.machValues = None
        self.totalTankVolume = None

    def setMass(self, mass):

        self.mass = mass

    def setDragCharacteristics(self, frontalArea, dragCurve):

        self.frontalArea = frontalArea

        # If the drag curve is given as a csv file name, it will be read into memory
        if type(dragCurve) is str:
            machList = []
            cdList = []

            file = openFromRelativePath(dragCurve)

            reader = csv.reader(file)
            for row in reader:
                machList.append(row[0])
                cdList.append(row[1])

            file.close()
            
            self.cdValues = machList
            self.machValues = cdList

            # Data is read in as strings, so this converts to floats
            for i in range(len(self.cdValues)):
                self.cdValues[i] = float(self.cdValues[i])
            for i in range(len(self.machValues)):
                self.machValues[i] = float(self.machValues[i])

        # If the drag curve is given as a list of two lists, they are assigned to instance variables
        elif type(dragCurve) is list:
            self.cdValues = float(dragCurve[0])
            self.machValues = float(dragCurve[1])
        else:
            sys.exit("dragCurve variable is not given in an accepted format")
    
    def setTotalTankVolume(self, totalTankVolume):

        self.totalTankVolume = totalTankVolume