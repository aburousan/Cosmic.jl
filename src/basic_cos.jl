using QuadGK




abstract type AbstractCosmology end
abstract type AbstractClosedCosmology <: AbstractCosmology end
abstract type AbstractFlatCosmology <: AbstractCosmology end
abstract type AbstractOpenCosmology <: AbstractCosmology end


# For flat ΛCDM
struct FlatLCDM{T <: Real} <: AbstractFlatCosmology
    h::T
    Ω_Λ::T
    Ω_m::T
    Ω_r::T
end
FlatLCDM(h::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real) = FlatLCDM(promote(float(h), float(Ω_Λ), float(Ω_m), float(Ω_r))...)
H_a2_by_H0(c::FlatLCDM, a) = sqrt(c.Ω_r + c.Ω_m * a + c.Ω_Λ * a^4)
# mostly we will use this

# For just an extension later will work on this
struct ClosedLCDM{T <: Real} <: AbstractClosedCosmology
    h::T
    Ω_k::T
    Ω_Λ::T
    Ω_m::T
    Ω_r::T
end
ClosedLCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real) = ClosedLCDM(promote(float(h), float(Ω_k), float(Ω_Λ), float(Ω_m),float(Ω_r))...)

struct OpenLCDM{T <: Real} <: AbstractOpenCosmology
    h::T
    Ω_k::T
    Ω_Λ::T
    Ω_m::T
    Ω_r::T
end
OpenLCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real) = OpenLCDM(promote(float(h), float(Ω_k), float(Ω_Λ), float(Ω_m),float(Ω_r))...)

function H_a2_by_H0(c::Union{ClosedLCDM,OpenLCDM}, a)
    a2 = a * a
    return sqrt(c.Ω_r + c.Ω_m * a + (c.Ω_k + c.Ω_Λ * a2) * a2)
end
# end here.

#For understanding w0 and wa see the paper.
#https://doi.org/10.1016/j.aop.2023.169244
for c in ("Flat", "Open", "Closed")
    name = Symbol("$(c)WCDM")
    @eval begin
        struct $(name){T <: Real} <: $(Symbol("Abstract$(c)Cosmology"))
            h::T
            Ω_k::T
            Ω_Λ::T
            Ω_m::T
            Ω_r::T
            w0::T
            wa::T
        end
        function $(name)(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real,
                         w0::Real, wa::Real)
            $(name)(promote(float(h), float(Ω_k), float(Ω_Λ), float(Ω_m),
                            float(Ω_r), float(w0), float(wa))...)
        end
    end
end
#For a more appropriate modle
function WCDM(h::Real, Ω_k::Real, Ω_Λ::Real, Ω_m::Real, Ω_r::Real, w0::Real, wa::Real)
    if Ω_k < 0
        ClosedWCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_r, w0, wa)
    elseif Ω_k > 0
        OpenWCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_r, w0, wa)
    else
        FlatWCDM(h, Ω_k, Ω_Λ, Ω_m, Ω_r, w0, wa)
    end
end

function H_a2_by_H0(c::Union{FlatWCDM,ClosedWCDM,OpenWCDM}, a)
    ade = exp((1 - 3 * (c.w0 + c.wa)) * log(a) + 3 * c.wa * (a - 1))
    return sqrt(c.Ω_r + (c.Ω_m + c.Ω_k * a) * a + c.Ω_Λ * ade)
end

# Final Cosmology Function
# All predefined values are taken from Baumann's book table-2.1
function cosmology(;h = 0.6774, Neff = 3.04, Ωk = 0, Ωm = 0.3089,
    Ωr = nothing, Tcmb = 2.7255, w0 = -1,wa = 0)
    if Ωr === nothing
        Ωγ = 4.48131e-7 * Tcmb^4 / h^2
        Ων = Neff * Ωγ * (7 / 8) * (4 / 11)^(4 / 3)
        Ωr = Ωγ + Ων
    end

    ΩΛ = 1 - Ωk - Ωm - Ωr

    if !(w0 == -1 && wa == 0)
        return WCDM(h, Ωk, ΩΛ, Ωm, Ωr, w0, wa)
    end

    if Ωk < 0
        return ClosedLCDM(h, Ωk, ΩΛ, Ωm, Ωr)
    elseif Ωk > 0
        return OpenLCDM(h, Ωk, ΩΛ, Ωm, Ωr)
    else
        return FlatLCDM(h, ΩΛ, Ωm, Ωr)
    end
end
#--------------------------------------------------------------------------------
#Conversion
Mpc_to_km(mp_km) = mp_km * 3.262e6 * 9.461e12
sec_to_year(sec_ye) = 3.171e-8*sec_ye
#--------------------------------------------------------------------------------
# Scale Factor
a(z) = 1/(1+z)#Sacle factor as a function of z
redshift(a) = 1/a - 1# Inverse function
#--------------------------------------------------------------------------------
# Hubble Rate
H_by_H0(c::AbstractCosmology, z) = (scal = a(z); H_a2_by_H0(c, scal) / scal^2)
H(c::AbstractCosmology, z) = 100 * c.h * H_by_H0(c, z)# in km /s / Mpc
#--------------------------------------------------------------------------------
# Distance
χ0(c::AbstractCosmology) = 2997.92458 / c.h# in Mpc
hubble_distance(c::AbstractCosmology, z) = χ0(c) / H_by_H0(c, z)# 1/H -> hubble radius
χ(c::AbstractCosmology, z::Real, ::Nothing; kws...) = QuadGK.quadgk(a->1 / H_a2_by_H0(c, a), a(z), 1; kws...)[1]
χ(c::AbstractCosmology, z₁::Real, z₂::Real; kws...) = QuadGK.quadgk(a->1 / H_a2_by_H0(c, a), a(z₂), a(z₁); kws...)[1]
# Times
η0(c::AbstractCosmology) = 3.262*9.46*1e18 / (1e2*c.h)# sec
hubble_time(c::AbstractCosmology, z) = η0(c) / H_by_H0(c, z)
T(c::AbstractCosmology, a0, a1; kws...) = QuadGK.quadgk(x->x / H_a2_by_H0(c, x), a0, a1; kws...)[1]
age(c::AbstractCosmology, z; kws...) = η0(c) * T(c, 0, a(z); kws...)
lookback_time(c::AbstractCosmology, z; kws...) = η0(c) * T(c, a(z), 1; kws...)
# Constants all are in MeV and C = 1, ħ = 1, kᵦ = 1
# H0 = (67.74/(3.26*9.5))*1e-18
# Ωm = 0.32
# ΩΛ = 0.68
# Ωγ = 5.35e-5
# Ων = 3.64e-5
# Ωr = Ωγ + Ων
# zeq = 3395
# All constants are from Baumann's book



function inte_for_age(a)
    deno = √(Ωr/a^2 + Ωm/a + ΩΛ*a^2)
    return 1/deno
end


age(a) = QuadGK.quadgk(inte_for_age, 0, a)[1]/H0 # In sec
my_f(x,y) = 2x+7y

my_g(x) = x^2

