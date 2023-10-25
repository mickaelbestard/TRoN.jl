export InitFunction
export Constant, set_value!
export Sinusoidal, set_Amplitude!, set_Phase!

abstract type InitFunction end 
(self::InitFunction)(v::Vector{Float64}) = self.(v) 

Base.@kwdef mutable struct Constant <: InitFunction
    value :: Float64 = 0.66
end

function (self::Constant)(::Float64)
    return self.value
end

function set_value!(self::Constant, value)
    self.value = value
end

# -----------------------------------------

Base.@kwdef mutable struct Sinusoidal <: InitFunction
    Amplitude :: Float64 = 0.3
    Phase :: Float64 = 1
end

function (self::Sinusoidal)(x::Float64)
    return self.Amplitude * ( sin(self.Phase * x) + one(x) )
end

function set_Amplitude!(self::Sinusoidal, value)
    self.Amplitude = value
end

function set_Phase!(self::Sinusoidal, value)
    self.Phase = value
end