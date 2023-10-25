export Diagnostics
### faire un dict à la place ?
struct Diagnostics 
    ### history : values for each iterations
    losshistory   :: Vector{Float64} 
    costhistory   :: Vector{Float64} 
    staffhistory  :: Vector{Float64} 
    totvarhistory :: Vector{Float64}  
    residualhistory :: Vector{Float64} 
    optimalityhistory :: Vector{Float64} # iter ↦ norm(Λ)
    descentstephistory :: Vector{Float64}
    ### hyperparam history
    thetaStaff_history :: Vector{Float64}
    thetaBV_history :: Vector{Float64}
    ### snapshots : values at optimal iteration 
    u_best :: Vector{Vector{Float64}}
    ρ_best :: Vector{Vector{Float64}}
    ### snapshots : values at last iteration 
    u_last :: Vector{Vector{Float64}}
    ρ_last :: Vector{Vector{Float64}}

    function Diagnostics(model::Model)
        losshistory = Float64[]
        costhistory = Float64[]
        staffhistory = Float64[]
        totvarhistory = Float64[]
        residualhistory = Float64[]
        optimalityhistory = Float64[]
        descentstephistory = Float64[]

        thetaStaff_history = Float64[]
        thetaBV_history = Float64[]

        u_best = 0. .* copy(model.u_t)
        ρ_best = 0. .* copy(model.ρ_t)

        u_last = 0. .* copy(model.u_t)
        ρ_last = 0. .* copy(model.ρ_t)

        new(losshistory, costhistory, staffhistory, totvarhistory, 
        residualhistory, optimalityhistory, descentstephistory,
        thetaStaff_history, thetaBV_history,
        u_best, ρ_best, u_last, ρ_last)
    end

    Diagnostics(a1::Vector{Float64}, a2::Vector{Float64}, 
                a3::Vector{Float64}, a4::Vector{Float64}, 
                a5::Vector{Float64}, a6::Vector{Float64}, 
                a7::Vector{Float64}, a8::Vector{Float64}, 
                a9::Vector{Float64}, 
                a10::Vector{Vector{Float64}}, 
                a11::Vector{Vector{Float64}},
                a12::Vector{Vector{Float64}}, 
                a13::Vector{Vector{Float64}}
                ) = new(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13)
end