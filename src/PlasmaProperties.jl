module PlasmaProperties

using Unitful
import PhysicalConstants.CODATA2018: c_0, ε_0, μ_0, k_B, m_e, e, m_u

export f_ce, f_pe, λ_D # Calculator functions
export eV, _m3, T, G, Hz, MHz # Unit shortcuts
export ε_0, μ_0, m_e, e, m_u # Physical constants, reexported from PhysicalConstants

const eV = u"eV"
const _m3 = u"m^-3"
const _cm3 = u"cm^-3"
const T = u"T"
const mT = u"mT"
const G = u"Gauss"
const Hz = u"Hz"
const MHz = u"MHz"
 
"""
    f_ce(B)

Compute the cyclotron frequency of electrons in MHz.
If B does not have Unitful units, it is assumed to be in Tesla.
"""
f_ce(B) = (e*B*1T/m_e /(2π) |> MHz)
f_ce(B::Unitful.BField) = (e*B/m_e /(2π) |> MHz)

"""
    f_pe(ne)

Compute the plasma frequency of electrons in MHz.
If ne does not have Unitful units, it is assumed to be in m^-3.
""" 
f_pe(ne) = (sqrt(ne*1u"m^-3"*e^2/(m_e*ε_0)) /(2π) |> MHz)
f_pe(ne::Quantity) = (sqrt(ne*e^2/(m_e*ε_0)) /(2π) |> MHz) 

"""
    λ_D(ne,Te)

Compute the Debye length of the plasm in m.
If ne does not have Unitful units, it is assumed to be in m^-3.
If Te does not have Unitful units, it is assumed to be in eV.
"""
λ_D(ne,Te) = sqrt(ε_0*Te*1u"eV"/(ne*1u"m^-3"*e^2)) |> u"m" 
λ_D(ne::Quantity,Te) = sqrt(ε_0*Te*1u"eV"/(ne*e^2)) |> u"m" 
λ_D(ne,Te::Unitful.Energy) = sqrt(ε_0*Te/(ne*1u"m^-3"*e^2)) |> u"m" 
λ_D(ne,Te::Unitful.Temperature) = sqrt(ε_0*k_B*Te/(ne*1u"m^-3"*e^2)) |> u"m" 
λ_D(ne::Quantity,Te::Unitful.Energy) = sqrt(ε_0*Te/(ne*e^2)) |> u"m" 
λ_D(ne::Quantity,Te::Unitful.Temperature) = sqrt(ε_0*k_B*Te/(ne*e^2)) |> u"m" 

end # module
