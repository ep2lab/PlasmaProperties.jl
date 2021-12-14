"""
PlasmaProperties, a simple plasma calculator module in Julia

The module contains basic functions to compute common plasma properties, such as characteristic frequencies and lengths.

It builds on Unitful, which provides a simple unit system for physical quantities.
"""
module PlasmaProperties

using Unitful 
import Unitful: Length, Time, Mass, Energy, Temperature, BField
using PhysicalConstants
me = PhysicalConstants.CODATA2018.m_e |> u"kg"
ε0 = PhysicalConstants.CODATA2018.ε_0 |> u"F/m"
μ0 = PhysicalConstants.CODATA2018.μ_0 |> u"N/A^2"
kB = PhysicalConstants.CODATA2018.k_B |> u"J/K"
qe = PhysicalConstants.CODATA2018.e |> u"C"
mu = PhysicalConstants.CODATA2018.m_u |> u"kg"

export Plasma # Plasma object type and constructors
export atomic_masses # Atomic masses
export f_ce, f_pe, rLe, λ_De # Calculator functions 
export ε0, μ0, me, qe # Physical constants, reexported from PhysicalConstants

@derived_dimension ParticleDensity Unitful.𝐋^-3
 
"""
check_or_assume(value,units)

If `value` is not a Unitful quantity, assume it is in the given `units`.
"""
function check_or_assume(value,units)
    if isa(value,Unitful.Quantity)
        value 
    else
        Unitful.Quantity(value,units)
    end
end

"""
    force_energy_units!(value)

If `value` is in temperature units, convert it to energy units.
"""
function force_energy_units!(value)
    if isa(value,Unitful.Temperature)
        value = kB*value 
    end
    value
end 

const atomic_masses = Dict([
        ("Ar", 39.948*mu),
        ("Xe", 131.293*mu),
        ])

const atomic_numbers = Dict([
        ("Ar", 1),
        ("Xe", 1),
        ])

"""
    P = Plasma()

Mutable type holding plasma properties like density, temperature, etc.
The show() method prints a table of the properties and derived quantities.

A name can be specified to quickly select the ion mass and charge. 
If units are not given, the default is `mi` in kg, `n` in 1/m^3, `Te` in eV, `B` in Gauss.
"""
mutable struct Plasma 
    name::String
    mi::Mass 
    Z::Int
    n::ParticleDensity 
    Te::Energy
    B::BField
end 
Plasma(mi,Z,n,Te,B) = Plasma("Custom",
                           check_or_assume(mi,"kg"),
                           Z,
                           check_or_assume(n,u"m^-3"),
                           force_energy_units!(check_or_assume(Te,u"eV")),
                           check_or_assume(B,u"Gauss"))
Plasma(name::String,n,Te,B) = Plasma(name,
                                     atomic_masses[name],
                                     atomic_numbers[name],
                                     check_or_assume(n,u"m^-3"),
                                     check_or_assume(Te,u"eV"),
                                     check_or_assume(B,u"Gauss"))

function Base.show(io::IO,P::Plasma)
    println(io,P.name," plasma with properties:")    
    println(io,"mi = ",P.mi)
    println(io,"Z = ",P.Z)
    println(io,"n = ",P.n)
    println(io,"Te = ",P.Te)
    println(io,"B = ",P.B)
    println(io,"f_ce = ",f_ce(P))
    println(io,"f_pe = ",f_pe(P))
    println(io,"rLe = ",rLe(P))
    println(io,"λ_De = ",λ_De(P))
end

f_c(m,Z,B) = Z*qe*check_or_assume(B,u"Gauss")/check_or_assume(m,u"kg") /(2π) |> u"MHz"

"""
    f_ce(B)
    f_ce(P::Plasma)

Compute the (unsigned) cyclotron frequency of electrons in MHz.
If `B` does not have Unitful units, it is assumed to be in Tesla.
"""
f_ce(B) = f_c(me,1,B)
f_ce(P::Plasma) = f_ce(P.B)

"""
    f_ci(mi,Z,B)
    f_ci(P::Plasma)

Compute the cyclotron frequency of ions of mass `mi` in MHz.
If `B` does not have Unitful units, it is assumed to be in Tesla.
"""
f_ci(mi,Z,B) = f_c(mi,Z,B)
f_ci(P::Plasma) = f_ci(P.mi,P.Z,P.B)

f_p(m,Z,n) = sqrt(check_or_assume(n,u"m^-3")*Z^2*qe^2/(check_or_assume(m,u"kg")*ε0)) /(2π) |> u"MHz"

"""
    f_pe(ne)
    f_pe(P::Plasma)

Compute the (unsigned) plasma frequency of electrons in MHz.
If `ne` does not have Unitful units, it is assumed to be in m^-3.
""" 
f_pe(ne) = f_p(me,1,ne)
f_pe(P::Plasma) = f_pe(P.n)

"""
    f_pi(mi,Z,ni)
    f_pi(P::Plasma)

Compute the plasma frequency of ions of mass `mi` in MHz.
If `ni` does not have Unitful units, it is assumed to be in m^-3.
""" 
f_pi(mi,Z,ni) = f_p(mi,Z,ni)
f_pi(P::Plasma) = f_pi(P.mi,P.Z,P.n)

rL(m,Z,B,vperp) = check_or_assume(m,u"kg")*check_or_assume(vperp,u"m/s")/(Z*qe*check_or_assume(B,u"Gauss")) |> u"m"

"""

    rLe(B,vperpe)
    
Compute the electron gyro radius in m.
If `B` does not have Unitful units, it is assumed to be in Tesla.
If `vperpe` does not have Unitful units, it is assumed to be in m/s.
"""
rLe(B,vperpe) = rL(me,1,B,vperpe)
rLe(P::Plasma) = rLe(P.B,sqrt(P.Te/me))

"""
    λ_De(ne,Te)
    λ_De(P::Plasma)

Compute the Debye length of the plasm in m.
If ne does not have Unitful units, it is assumed to be in m^-3.
If Te does not have Unitful units, it is assumed to be in eV.
"""
λ_De(ne,Te) = sqrt(ε0*force_energy_units!(check_or_assume(Te,u"eV"))/(check_or_assume(ne,u"m^-3")*qe^2)) |> u"m"  
λ_De(P::Plasma) = λ_De(P.n,P.Te)

end # module
