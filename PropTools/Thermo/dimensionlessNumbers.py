def reynoldsNumber(density: float, velocity: float, characteristicLength: float, visocity: float) -> float:

    return (density * velocity * characteristicLength) / visocity

def prandtlNumber(specificHeat: float, viscosity: float, thermalConductivity: float) -> float:

    return (specificHeat * viscosity) / thermalConductivity