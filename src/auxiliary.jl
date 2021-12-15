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