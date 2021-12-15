"""
PlasmaProperties, a simple plasma calculator module in Julia

The module contains basic functions to compute common plasma properties, such as characteristic frequencies and lengths.

It builds on Unitful, which provides a simple unit system for physical quantities.
"""
module PlasmaProperties

using Unitful 
import Unitful: Length, Time, Mass, Energy, Temperature, BField, EField
using PhysicalConstants
ε0 = PhysicalConstants.CODATA2018.ε_0 |> u"F/m"
μ0 = PhysicalConstants.CODATA2018.μ_0 |> u"N/A^2"
me = PhysicalConstants.CODATA2018.m_e |> u"kg"
qe = PhysicalConstants.CODATA2018.e |> u"C"
kB = PhysicalConstants.CODATA2018.k_B |> u"J/K"
mu = PhysicalConstants.CODATA2018.m_u |> u"kg"

export QuasineutralPlasma # Plasma object type and constructors
export atomic_masses # Atomic masses
export f_ce, f_pe, r_Le, r_Lsi, λ_De, c_s, c_the, ν_ei, ν_ie, ν_ee, χ # Calculator functions 
export ε0, μ0, me, qe, kB # Physical constants, reexported from PhysicalConstants
 
include("auxiliary.jl")
include("types.jl")
include("functions.jl")

end # module
