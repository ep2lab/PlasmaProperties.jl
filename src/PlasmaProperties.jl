module PlasmaProperties

import PhysicalConstants.CODATA2018: c_0, ε_0, μ_0, k_B, m_e, e, m_u
using Unitful

export f_ce, f_pe

"""
    f_ce(B)

    Compute the cyclotron frequency of electrons in MHz.
"""
f_ce(B) = (e*B/m_e /(2π) |> u"MHz")


"""
    f_pe(B)

    Compute the plasma frequency of electrons in MHz.
""" 
f_pe(ne) = (sqrt(ne*e^2/(m_e*ε_0)) /(2π) |> u"MHz")


end # module
