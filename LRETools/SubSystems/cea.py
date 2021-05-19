from rocketcea.cea_obj_w_units import CEA_Obj
import sys

# Iterates through mixture ratios, performing CEA on each one, to find the ratio that produces the highest specific impulse

def idealMixtureRatioCEA(fuel, ox, chamberPressure, ambientPressure=1.01325, condition='equilibrium', throatCondition='equilibrium', searchResolution=0.1, startSearch=2):

    CEA = CEA_Obj(oxName=ox, fuelName=fuel, pressure_units='bar')

    # Converts equlibirum/frozen statements into a flag for the CEA functions
    if condition == 'equilibrium':
        c = 0
    elif condition == 'frozen':
        c = 1
    else:
        sys.exit('CEA condition can only be "equilibrium" or "frozen"')
    
    if throatCondition == 'equilibrium':
        ct = 0
    elif throatCondition == 'frozen':
        ct = 1
    else:
        sys.exit('CEA throat condition can only be "equilibrium" or "frozen"')

    # Initialises variables for following loop
    i = startSearch
    I_sp = 1
    I_sp_prev = 0
    expansionRatio = 0

    # Performs mixture ratio iterations until the specific impulse starts dropping
    while I_sp >= I_sp_prev:

        I_sp_prev = I_sp

        expansionRatio = CEA.get_eps_at_PcOvPe(Pc=chamberPressure, MR=i, PcOvPe=chamberPressure/ambientPressure, frozen=c, frozenAtThroat=ct)
        I_sp = CEA.get_Isp(Pc=chamberPressure, MR=i, eps=expansionRatio, frozen=c, frozenAtThroat=ct)
        chamberTemp = CEA.get_Tcomb(Pc=chamberPressure, MR=i)

        i = i + searchResolution
    
    # Rewinds back one iteration to get the ideal specific impulse mixture ratio
    idealMixtureRatio = i - searchResolution
    expansionRatio = CEA.get_eps_at_PcOvPe(Pc=chamberPressure, MR=idealMixtureRatio, PcOvPe=chamberPressure/ambientPressure, frozen=c, frozenAtThroat=ct)
    I_sp = CEA.get_Isp(Pc=chamberPressure, MR=idealMixtureRatio, eps=expansionRatio, frozen=c, frozenAtThroat=ct)
    chamberTemp = CEA.get_Tcomb(Pc=chamberPressure, MR=i)
    speciesMassFractions = CEA.get_SpeciesMassFractions(Pc=chamberPressure, MR=idealMixtureRatio, eps=expansionRatio, frozen=c, frozenAtThroat=ct, min_fraction=0.005)
    throatMolWtGamma = CEA.get_Throat_MolWt_gamma(Pc=chamberPressure, MR=idealMixtureRatio, eps=expansionRatio, frozen=ct)

    # Outputs a dictionary containing useful parameters of the ideal mixture ratio performance
    output = {
        "Mixture Ratio": idealMixtureRatio,
        "Specific Impulse": I_sp,
        "Chamber Temperature": chamberPressure,
        "Expansion Ratio": expansionRatio,
        "Species Mass Fractions": speciesMassFractions,
        "Throat Molecular Weight": throatMolWtGamma[0],
        "Gamma": throatMolWtGamma[1]
    }

    return(output)
