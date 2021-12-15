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