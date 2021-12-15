# PlasmaProperties.jl

A simple calculator of plasma properties like characteristic lengths and frequencies.


[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://ep2lab.github.io/PlasmaProperties.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://ep2lab.github.io/PlasmaProperties.jl/dev)
[![Build status (Github Actions)](https://github.com/ep2lab/PlasmaProperties.jl/workflows/CI/badge.svg)](https://github.com/ep2lab/PlasmaProperties.jl/actions)

## Usage

Load the `PlasmaProperties` package and the `Unitful` package:

```julia
using PlasmaProperties, Unitful
```

Create a QuasineutralPlasma object. Use a name 'Xe' or 'Ar' to set the ion mass to the correct value for xenon and argon, respectively. Singly-charged ions are assumed. Other values of mass and charge state can be easily added to the package. You can give the values of the plasma density, neutral density, electron temperature, and the fields at creation time, or change them afterwards:

```julia
P = QuasineutralPlasma('Xe', 1e18u"m^-3", 3u"eV")
P.B = 200u"Gauss"
```

The units at creation time are optional; the usual ones are assumed if none are given.

The overloaded `Base.show()` method displays all the information of the plasma when printing `P` in the REPL. You can also use functions to compute a given quantity. E.g. for the electron cyclotron frequency in GHz,

```julia
f_ce(P)
```
