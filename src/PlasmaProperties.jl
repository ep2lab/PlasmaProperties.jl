module PlasmaProperties

import PhysicalConstants.CODATA2018: c_0, ε_0, μ_0, k_B, m_e, e, m_u
using Unitful

"""
    ω_ce(B)

    Compute the cyclotron frequency of electrons in units of Hz.
"""
ω_ce(B) = (e*B/(2π*m_e) |> u"Hz")

end # module
