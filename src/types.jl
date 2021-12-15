const atomic_masses = Dict([
        ("Ar", 39.948*mu),
        ("Xe", 131.293*mu),
        ])

const atomic_numbers = Dict([
        ("Ar", 1),
        ("Xe", 1),
        ])


@derived_dimension ParticleDensity Unitful.𝐋^-3

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