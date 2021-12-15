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
The plasma is quasineutral, cold-ion.
The show() method prints a table of the properties and derived quantities.

A name can be specified to quickly select the ion mass and charge. 
If units are not given, the default is `mi` in kg, `n` in 1/m^3, `Te` in eV, `B` in Gauss.
"""
mutable struct QuasineutralPlasma 
    name::String

    np::ParticleDensity 
    Te::Energy

    B::BField
    E::EField
    nn::ParticleDensity

    mi::Mass 
    Zi::Int
end  
QuasineutralPlasma(name::String="Xe",np=1e18u"m^-3",Te=5u"eV";B=0u"Gauss",E=0u"V/m",nn=0u"m^-3") = 
    QuasineutralPlasma(
    name,
    check_or_assume(np,u"m^-3"),
    check_or_assume(Te,u"eV"),
    check_or_assume(B,u"Gauss"),
    check_or_assume(E,u"V/m"),
    check_or_assume(nn,u"m^-3"),
    atomic_masses[name],
    atomic_numbers[name],
    )

function Base.show(io::IO,P::QuasineutralPlasma)
    println(io,P.name," quasineutral plasma")
    println(io,"")

    println(io,"Plasma density:")    
    println(io,"n = ",P.np)
    println(io,"")

    println(io,"Field properties:")
    println(io,"B = ",P.B)
    println(io,"E = ",P.E)
    println(io,"")

    println(io,"Cold ion properties:")
    println(io,"mi = ",P.mi)
    println(io,"Zi = ",P.Zi," (charge number)")
    println(io,"c_s = ",c_s(P)," (sound speed based on electron temperature)")
    println(io,"f_ci = ",f_ci(P))
    println(io,"r_Lsi = ",r_Lsi(P)," (ion gyroradius based on sound speed)")
    println(io,"")

    println(io,"Electron properties:")
    println(io,"Te = ",P.Te)
    println(io,"c_the = ",c_the(P))
    println(io,"f_ce = ",f_ce(P))
    println(io,"r_Le = ",r_Le(P))
    println(io,"f_pe = ",f_pe(P))
    println(io,"λ_De = ",λ_De(P))
    println(io,"")

    println(io,"Cold neutral properties:")
    println(io,"nn = ",P.nn)
    println(io,"")

    println(io,"Collisional properties:")
    println(io,"ν_ei = ",ν_ei(P))
    println(io,"ν_ie = ",ν_ie(P))
    println(io,"ν_ee = ",ν_ee(P))
    println(io,"χ = ",χ(P)," (Hall parameter based on ν_ei only)")
    println(io,"")
end                                     