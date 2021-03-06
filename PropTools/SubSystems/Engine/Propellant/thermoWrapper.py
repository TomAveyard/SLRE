import sys
from thermo import ChemicalConstantsPackage, PRMIX, CEOSLiquid, CEOSGas, FlashPureVLS, FlashVL
import thermo
from thermo.interaction_parameters import IPDB

class thermoWrapper:

    def __init__(self, name, massFractions=None):
        
        self.name = name
        self.massFractions = massFractions

        # Initialise variables that will be assigned later
        self.constants = None
        self.properties = None
        self.result = None

        # Loads a pure substance if a single string is given
        if type(name) is str:
            self.type = "Pure"
            self.flasher = self.loadPureSubstance()
            self.moleFractions = None

            self.averageMW = self.constants.MWs[0]

        # Loads a mixture if a list of strings is given
        elif type(name) is list:
            self.type = "Mixture"
            self.flasher = self.loadMixture()

            if sum(self.massFractions) != 1:
                sys.exit("Invalid: Mass fractions do not add up to 1")
            
            self.moleFractions = self.massToMoleFractions()

            # Calculates average molecular weight
            self.averageMW = 0
            for i in range(len(self.constants.MWs)):

                self.averageMW = self.constants.MWs[i] * self.moleFractions[i]

    # If a pure substance is given, then the thermo flasher object for pure substances is created and stored
    def loadPureSubstance(self):

        # Load the constant properties and correlation properties
        self.constants, self.properties = ChemicalConstantsPackage.from_IDs([self.name])

        # Configure the liquid and gas phase objects
        eos_kwargs = dict(Tcs=self.constants.Tcs, Pcs=self.constants.Pcs, omegas=self.constants.omegas)
        liquid = CEOSLiquid(PRMIX, HeatCapacityGases=self.properties.HeatCapacityGases, eos_kwargs=eos_kwargs)
        gas = CEOSGas(PRMIX, HeatCapacityGases=self.properties.HeatCapacityGases, eos_kwargs=eos_kwargs)

        # Create a flash object with possible phases of 1 gas and 1 liquid
        flasher = FlashPureVLS(self.constants, self.properties, gas=gas, liquids=[liquid], solids=[])

        return(flasher)

    # If a mixture is given, then the thermo flasher object for mixtures is created and stored
    def loadMixture(self):
        
        # Load the constant properties and correlation properties
        self.constants, self.properties = ChemicalConstantsPackage.from_IDs(self.name)

        # Configure the liquid and gas phase objects
        kijs = IPDB.get_ip_asymmetric_matrix('ChemSep PR', self.constants.CASs, 'kij')
        eos_kwargs = {'Pcs': self.constants.Pcs, 'Tcs': self.constants.Tcs, 'omegas': self.constants.omegas, 'kijs': kijs}
        gas = CEOSGas(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=self.properties.HeatCapacityGases)
        liquid = CEOSLiquid(PRMIX, eos_kwargs=eos_kwargs, HeatCapacityGases=self.properties.HeatCapacityGases)

        # Create a flash object with possible phases of 1 gas and 1 liquid
        flasher = FlashVL(self.constants, self.properties, liquid=liquid, gas=gas)

        return(flasher)

    # Converts the given mass fractions to mole fractions to give to thermo
    def massToMoleFractions(self):

        molecularWeights = self.constants.MWs
        moles = []

        for i in range(len(self.massFractions)):

            moles.append(self.massFractions[i] / molecularWeights[i])

        totalMoles = sum(moles)

        for i in range(len(moles)):

            moles[i] = moles[i] / totalMoles

        return(moles)

    def densityToMolarVolume(self, D):

        molarVolume = self.averageMW / (D * 1000)
        return(molarVolume)

    def molarVolumeToDensity(self, molarVolume):

        D = (self.averageMW) / (molarVolume * 1000)

        return(D)

    def perMassToPerMol(self, quantityPerMol):

        return quantityPerMol * self.averageMW

    def calculate(self, T=None, P=None, S=None, H=None, D=None, Q=None):
        
        # If density is given, it is converted to molar volume as that is what needs to be given to thermo
        if D != None:
            molarVolume = self.densityToMolarVolume(D)
        else:
            molarVolume = None

        # If entropy or enthalpy are given, they are converted to per mol rather than per kg, as that is what thermo requires
        if S != None:
            S = self.perMassToPerMol(S) / 1000
        if H != None:
            H = self.perMassToPerMol(H) / 1000

        # Calls the thermo flash function and stores results
        if self.type == "Pure":

            #print(T, P, S, H, molarVolume, Q)
            result = self.flasher.flash(T=T, P=P, S=S, H=H, V=molarVolume, VF=Q)

        elif self.type == "Mixture":
            
            result = self.flasher.flash(T=T, P=P, S=S, H=H, V=molarVolume, VF=Q, zs=self.moleFractions)

        self.T = result.T
        self.P = result.P
        self.S = result.S_mass()
        self.H = result.H_mass()
        self.D = self.molarVolumeToDensity(result.V())
        self.Q = result.VF
        self.Cp = result.Cp_mass()
        self.Cv = result.Cv_mass()
        self.gamma = self.Cp/self.Cv
        self.mu = result.mu()
        self.kappa = result.k()