<div align="center">
    <img src="Images/slreLogo.png" alt="SLRE Logo" width="400"/>
</div>
<div align="center">
    <p><b>S</b>imulation of <b>L</b>iquid <b>R</b>ocket <b>E</b>ngines</p>
    <p><i>(Pronounced "Celery")</i></p>
</div>

---

SLRE is a Python package that provides simulation tools for the preliminary design of liquid rocket engines.

---

## Features

Features include:
- Integrated CEA analysis through [ceapy.py](https://docs.pycea.org/chapters/introduction.html)
- Integrated thermodynamic property querying through [CoolProp](http://www.coolprop.org/)
- Generation of nozzle geometry (conical or Rao bell nozzle approximation)
- Generation of combustion chamber geometry
- One-dimensional regenerative cooling analysis
- Engine cycle analysis

---

## Installation

I have largely halted development of SLRE, but may revisit it in the future. 
The code is functional, but I have not extensively tested for bugs and usability, so be aware of this if you decide to use it.

To install:

1.  Clone the repository:

```shell
git clone https://github.com/TomAveyard/SLRE.git
```

2. Create a python virtual environment. If using VS Code, follow the following instructions: https://code.visualstudio.com/docs/python/environments

3. Install the requirements using:

```shell
pip install -r requirements.txt
```

If using VS Code you may have already been prompted to do this when setting up the virtual environment

4. You can test that everything is running correctly by running testScript.py in the main folder. This should result in an engine being solved for and the following output:

```
----------
Engine Performance
----------
Specific Impulse: 217.25 s
Exit Velocity: 2131.23 m/s
Thrust Coefficient: 1.39
C Star: 1552.26 m/s
----------
Mass Flow Rates
----------
Fuel Mass Flow Rate: 1.04 kg/s
Oxidiser Mass Flow Rate: 3.65 kg/s
Total Mass Flow Rate: 4.69 kg/s
----------
Temperatures
----------
Chamber Temperature: 2812.47 K
Throat Temperature: 2545.13 K
Exit Temperature: 1579.48 K
----------
Sizes
----------
Chamber Radius: 68.1 mm
Chamber Area: 14570.48 mm^2
Throat Radius: 34.05 mm
Throat Area: 3642.62 mm^2
Exit Radius: 62.11 mm
Exit Area: 12118.89 mm^2
----------
Power Balance
----------
Regenerative Cooling Heat Power: 1267.47 kW
Fuel Pump Power: -9.86 kW
Total Pump Power: -9.86 kW
Turbine Power: 10.21 kW
Power Balance: 0.36 kW
----------
Fuel Line Component Outputs
----------
----------
Fuel Pipe 1 Results
----------
Substance: propanol
Mass Flow Rate: 1.04 kg/s
Inlet Temperature: 298.0 K
Inlet Pressure: 3.0 Bar
Outlet Temperature: 298.0 K
Outlet Pressure: 2.91 Bar
Pressure Loss: 0.09 Bar
----------
----------
Fuel Pump Results
----------
Substance: propanol
Mass Flow Rate: 1.04 kg/s
Inlet Temperature: 298.0 K
Inlet Pressure: 2.91 Bar
Outlet Temperature: 300.27 K
Outlet Pressure: 40.0 Bar
Temperature Rise: -2.27 K
Pressure Rise: 37.09 Bar
Head Rise: 482.3 m
Power: -9.86 kW
----------
----------
Fuel Pipe 2 Results
----------
Substance: propanol
Mass Flow Rate: 1.04 kg/s
Inlet Temperature: 300.27 K
Inlet Pressure: 40.0 Bar
Outlet Temperature: 300.27 K
Outlet Pressure: 39.93 Bar
Pressure Loss: 0.07 Bar
----------
----------
Fuel Regen Cooling Results
----------
Substance: propanol
Mass Flow Rate: 1.04 kg/s
Inlet Temperature: 300.27 K
Inlet Pressure: 39.93 Bar
Outlet Temperature: 581.9 K
Outlet Pressure: 39.28 Bar
Temperature Rise: 281.64 K
Pressure Loss: 0.65 Bar
Max Heat Flux: 8103.84 kW/m^2
Max Gas Side Wall Temp: 660.98 K
----------
----------
Fuel Pipe 3 Results
----------
Substance: propanol
Mass Flow Rate: 1.04 kg/s
Inlet Temperature: 581.9 K
Inlet Pressure: 39.28 Bar
Outlet Temperature: 581.9 K
Outlet Pressure: 35.93 Bar
Pressure Loss: 3.35 Bar
----------
----------
Fuel Turbine Results
----------
Substance: propanol
Mass Flow Rate: 1.04 kg/s
Inlet Temperature: 581.9 K
Inlet Pressure: 35.93 Bar
Outlet Temperature: 570.05 K
Outlet Pressure: 26.51 Bar
Temperature Loss: 11.86 K
Pressure Loss: 9.42 Bar
Power: 10.21 kW
----------
----------
Fuel Pipe 4 Results
----------
Substance: propanol
Mass Flow Rate: 1.04 kg/s
Inlet Temperature: 570.05 K
Inlet Pressure: 26.51 Bar
Outlet Temperature: 570.05 K
Outlet Pressure: 21.26 Bar
Pressure Loss: 5.25 Bar
----------
----------
Oxidiser Line Component Outputs
----------
----------
Oxidiser Pipe 1 Results
----------
Substance: nitrous oxide
Mass Flow Rate: 3.65 kg/s
Inlet Temperature: 253.0 K
Inlet Pressure: 25.0 Bar
Outlet Temperature: 253.0 K
Outlet Pressure: 24.87 Bar
Pressure Loss: 0.13 Bar
----------
```

5. Any custom scripts using the package should be placed in and run from the main folder. testScript.py should serve as a useful guide for constructing your scripts.

---

## Example of Usage

SLRE has been designed to be (hopefully) quick and easy to use.

To begin, a thrust chamber object is created:

```Python
engineThrustChamber = ThrustChamber(fuelName='methane', oxName='oxygen', thrust=10*10**3, chamberPressure=40, fac=True, contractionRatio=3, ambientPressure=1.01325)
```

The combustion chamber and nozzle geometry for the thrust chamber is then generated, then combined into the overall thrust chamber geometry:

```Python
engineThrustChamber.getChamberGeometry(lStar=2.25,
                                       contractionLength=0.065, 
                                       numberOfPointsConverging=100,
                                       numberOfPointsStraight=50)
testThrustChamber.getRaoBellNozzleGeometry(0.8, numberOfPoints=100)
testThrustChamber.getThrustChamberCoords()
```

This engine involves regenerative cooling so parameters for the solver are now defined:

```Python
solverParameters = SolverParameters(bartzEquationCoefficient=0.026, coolantSideHeatTransferCorrelation="sieder-tate")
```

To perform a regenerative cooling analysis, the cooling channels are defined, then fed into a regenerative cooling object along with the thrust chamber and solver paramters:

```Python
engineCoolingChannels = CoolingChannels(numberOfChannels=96, wallThickness=0.9e-3, midRibThickness=1.5e-3, channelHeight=1.25e-3, wallConductivity=365, wallRoughnessHeight=1e-6)
engineRegenerativeCooling = RegenerativeCooling(thrustChamber=testThrustChamber, coolingChannels=testCoolingChannels, solverParameters=solverParameters)
```

To perform a cycle analysis, the state of the propellants in the tank are defined and Tank objects created to represent them:

```Python
fuel = Propellant(testThrustChamber.fuel.name)
ox = Propellant(testThrustChamber.ox.name)

fuel.defineState("T", 108, "P", 3*10**5)
ox.defineState("T", 60, "P", 3*10**5)

fuelTank = Tank(fuel)
oxTank = Tank(ox)
```

Other components can be defined e.g:

```Python
fuelPump = Pump(isentropicEfficiency=0.7, outletPressure=70e5)
turbine = Turbine(isentropicEfficiency=0.7, outletPressure=testThrustChamber.injectionPressure*10**5)

oxPump = Pump(isentropicEfficiency=0.7, outletPressure=testThrustChamber.injectionPressure*10**5)
```

These can then be combined together to form their respective lines in the cycle:

```Python
fuelLine = Line(inletState=fuelTank.outletState, massFlowRate=engineThrustChamber.fuelMassFlowRate, components=[fuelPump, engineRegenerativeCooling, turbine])
oxLine = Line(inletState=oxTank.outletState, massFlowRate=thrustChamber.oxMassFlowRate, components=[oxPump])
```

The cycle is then formed from these lines:

```Python
cycle = Cycle(fuelLine=fuelLine, oxLine=oxLine, thrustChamber=engineThrustChamber)
```

---
