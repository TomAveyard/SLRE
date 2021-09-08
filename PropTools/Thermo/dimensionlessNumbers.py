def reynoldsNumber(density, velocity, characteristicLength, visocity):

    return (density * velocity * characteristicLength) / visocity

def prandtlNumber(specificHeat, viscosity, thermalConductivity):

    return (specificHeat * viscosity) / thermalConductivity