"""
    f_c(m,Z,B) 

Compute the (unsigned) cyclotron frequency of a charged species in MHz.
If `m` does not have units, it is assumed to be in kg.
If `B` does not have Unitful units, it is assumed to be in Tesla.
"""
f_c(m,Z,B) = abs(Z*qe*check_or_assume(B,u"Gauss")/check_or_assume(m,u"kg") /(2π)) |> u"MHz"

"""
    f_pe(m,Z,n) 

Compute the plasma frequency of a charged species in MHz.
If `m` does not have units, it is assumed to be in kg.
If `n` does not have Unitful units, it is assumed to be in m^-3.
""" 
f_p(m,Z,n) = sqrt(check_or_assume(n,u"m^-3")*Z^2*qe^2/(check_or_assume(m,u"kg")*ε0)) /(2π) |> u"MHz"

"""
    rL(m,Z,B,vperp)
    
Compute the gyro radius of a charged species in m.
If `m` does not have units, it is assumed to be in kg.
If `B` does not have Unitful units, it is assumed to be in Tesla.
If `vperp` does not have Unitful units, it is assumed to be in m/s.
"""
r_L(m,Z,B,vperp) = abs(check_or_assume(m,u"kg")*check_or_assume(vperp,u"m/s")/(Z*qe*check_or_assume(B,u"Gauss"))) |> u"m"


"""
    cth(m,Z,T)
    
Compute the thermal speed of a charged species in m/s.
This ignores numerical factors and gammae.
"""
c_th(m,Z,T) = sqrt(Z*force_energy_units!(check_or_assume(T,u"eV"))/check_or_assume(m,u"kg")) |> u"m/s"

# ----------------------------------------------------------------------

"""
    f_ce(B)
    f_ce(P::QuasineutralPlasma)

Compute the (unsigned) cyclotron frequency of electrons in GHz.
If `B` does not have Unitful units, it is assumed to be in Tesla.
"""
f_ce(B) = f_c(me,1,B) |> u"GHz"
f_ce(P::QuasineutralPlasma) = f_ce(P.B)

"""
    f_ci(mi,Zi,B)
    f_ci(P::QuasineutralPlasma)

Compute the (unsigned) cyclotron frequency of ions of mass `mi` in kHz.
If `B` does not have Unitful units, it is assumed to be in Tesla.
"""
f_ci(mi,Zi,B) = f_c(mi,Zi,B) |> u"kHz"
f_ci(P::QuasineutralPlasma) = f_ci(P.mi,P.Zi,P.B)

"""
    f_pe(ne)
    f_pe(P::QuasineutralPlasma)

Compute the plasma frequency of electrons in GHz.
If `ne` does not have Unitful units, it is assumed to be in m^-3.
""" 
f_pe(ne) = f_p(me,1,ne) |> u"GHz"
f_pe(P::QuasineutralPlasma) = f_pe(P.np)

"""
    f_pi(mi,Z,ni)
    f_pi(P::QuasineutralPlasma)

Compute the plasma frequency of ions of mass `mi` in MHz.
If `ni` does not have Unitful units, it is assumed to be in m^-3.
""" 
f_pi(mi,Zi,ni) = f_p(mi,Zi,ni) |> u"kHz"
f_pi(P::QuasineutralPlasma) = f_pi(P.mi,P.Zi,P.np)

"""
    r_Le(B,vperpe)
    r_Le(P::QuasineutralPlasma)
    
Compute the electron gyro radius in mm.
If `B` does not have Unitful units, it is assumed to be in Tesla.
If `vperpe` does not have Unitful units, it is assumed to be in m/s.
"""
r_Le(B,vperpe) = r_L(me,1,B,vperpe) |> u"mm"
r_Le(P::QuasineutralPlasma) = r_Le(P.B,sqrt(P.Te/me))

"""
    r_Lsi(mi,Z,B,Te)
    r_Lsi(P::QuasineutralPlasma)
    
Compute the sonic ion gyro radius based on the sound speed, in cm.
If `mi` does not have units, it is assumed to be in kg.
If `B` does not have Unitful units, it is assumed to be in Tesla.
If `Te` does not have Unitful units, it is assumed to be in eV.
"""
r_Lsi(mi,Zi,B,Te) = r_L(mi,Zi,B,sqrt(force_energy_units!(check_or_assume(Te,u"eV"))/check_or_assume(mi,u"kg"))) |> u"cm"
r_Lsi(P::QuasineutralPlasma) = r_Lsi(P.mi,P.Zi,P.B,P.Te)

"""
    λ_De(ne,Te)
    λ_De(P::QuasineutralPlasma)

Compute the Debye length of the plasm in mm.
If ne does not have Unitful units, it is assumed to be in m^-3.
If Te does not have Unitful units, it is assumed to be in eV.
"""
λ_De(ne,Te) = sqrt(ε0*force_energy_units!(check_or_assume(Te,u"eV"))/(check_or_assume(ne,u"m^-3")*qe^2)) |> u"mm"  
λ_De(P::QuasineutralPlasma) = λ_De(P.np,P.Te)

"""
    c_s(mi,Zi,Te)
    c_s(P::QuasineutralPlasma)
    
Compute the sound speed of ions in m/s.
This ignores numerical factors and gammae.
"""
c_s(mi,Zi,Te) = c_th(mi,Zi,Te) # Just an alias for cth
c_s(P::QuasineutralPlasma) = c_s(P.mi,P.Zi,P.Te)

"""
    c_the(Te)
    c_the(P::QuasineutralPlasma)
    
Compute the thermal speed of electrons in km/s.
This ignores numerical factors and gammae.
"""
c_the(Te) = c_th(me,1,Te) |> u"km/s"
c_the(P::QuasineutralPlasma) = c_the(P.Te)


"""
    CoulombLog(ne,Te)

Approximate the Coulomb logarithm.
"""
CoulombLog(np,Te) = 9+0.5*log(1e18u"m^-3"/np * Te^3/1u"eV^3") # approximation  

"""
    ν_ei(Zi,np,Te)
    ν_ei(P::QuasineutralPlasma)

Electron-ion momentum collision frequency in MHz.
"""
function ν_ei(Zi,np,Te) 
    np = check_or_assume(np,u"m^-3")
    Te = force_energy_units!(check_or_assume(Te,u"eV"))
    return sqrt(2) * np * Zi^2 * qe^4 * CoulombLog(np,Te) /
        (12π^(3/2) * ε0^2 * sqrt(me) * Te^(3/2)) |> u"MHz"
end
ν_ei(P::QuasineutralPlasma) = ν_ei(P.Zi,P.np,P.Te)

"""
    ν_ie(mi,Zi,np,Te)
    ν_ie(P::QuasineutralPlasma)

Ion-electron momentum collision frequency in MHz.
"""
ν_ie(mi,Zi,np,Te) = me/mi*ν_ei(Zi,np,Te)
ν_ie(P::QuasineutralPlasma) = ν_ie(P.mi,P.Zi,P.np,P.Te)

"""
    ν_ee(ne,Te)
    ν_ee(P::QuasineutralPlasma)

Electron-electron momentum collision frequency in MHz.
"""
function ν_ee(ne,Te) 
    ne = check_or_assume(ne,u"m^-3")
    Te = force_energy_units!(check_or_assume(Te,u"eV"))    
    return ne * qe^4 * CoulombLog(ne,Te) /
        (12π^(3/2) * ε0^2 * sqrt(me) * Te^(3/2)) |> u"MHz"
end                 
ν_ee(P::QuasineutralPlasma) = ν_ee(P.np,P.Te)

"""
    χ(Zi,np,Te,B)
    χ(P::QuasineutralPlasma)

Hall parameter based on the electron-ion collision frequency only.
"""
χ(Zi,np,Te,B) = f_ce(B)/ν_ei(Zi,np,Te) |> Unitful.NoUnits
χ(P::QuasineutralPlasma) = χ(P.Zi,P.np,P.Te,P.B)