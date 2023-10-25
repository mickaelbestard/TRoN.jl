mutable struct Junction
    #TODOOO : ajouter coordonnées pour visualisation automatique ?
    Nin        :: Int64            # number of ingoing roads
    Nout       :: Int64            # number of outoing roads
    InNeigh    :: Vector{Int64}    # incoming edges
    OutNeigh   :: Vector{Int64}    # outgoing edges
    Jmat       :: Matrix{Union{Real, ForwardDiff.Dual}}   # transition matrix (dynamic : α(ub, uc))
    α          :: Matrix{Union{Real, ForwardDiff.Dual}}   # transition matrix (static : α∈R)
    priorities :: Vector{Float64}  # vector of priorities
    γin_max :: Vector{Real}
    γout_max :: Vector{Real}

    function Junction()
        Nin = 0
        Nout = 0
        InNeigh  = Int64[]
        OutNeigh = Int64[]
        Jmat     = zeros(Float64, (Nout, Nin))
        α        = zeros(Float64, (Nout, Nin))
        priorities = Float64[]
        γin_max = Vector{Real}(undef, 2)
        γout_max = Vector{Real}(undef, 2)
        new(Nin, Nout, InNeigh, OutNeigh, Jmat, α, priorities, γin_max, γout_max)
    end
end # class Junction


# ----- methods of class Junction -----

function set_priorities(self::Junction, v::Vector{Float64})
    self.priorities = v
end

function set_in_data(self::Junction, neigh::Int64)
    push!(self.InNeigh, neigh)
    self.Nin += 1
end

function set_out_data(self::Junction, neigh::Int64)
    push!(self.OutNeigh, neigh)
    self.Nout += 1
end

function set_Jmat(self::Junction, Jmat)
    self.Jmat = Jmat
    self.α    = copy(Jmat)
end

# call 1x1, 1x2, 2x1, 2x2 
function LP(self::Junction, γin_max::Vector{T}, γout_max::Vector{R}, uJ::Vector{S}) where {T,R,S <: Real}
# function LP(self::Junction, γin_max, γout_max, uJ) 
    # println("@LP()")

    ### create buffers
    γin  = copy(γin_max)
    γout = copy(γout_max)

    ### find the type of junction among {nxm; n,m ∈ {1,2}}
    if length(γin_max) == 1
        if length(γout_max) == 1 #? 1x1 

            LP_1x1!(γin, γin_max, γout_max, uJ)

        elseif length(γout_max) == 2 #? 1x2

            LP_1x2!(self, γin, γin_max, γout_max, uJ)

        end

    elseif length(γin_max) == 2

        if length(γout_max) == 1 #? 2x1, A=[1 1]

            LP_2x1!(self, γin, γin_max, γout_max, uJ)

        elseif length(γout_max) == 2 #? 2x2, A=[α β; (1-α) (1-β)] 

            LP_2x2!(self, γin, γin_max, γout_max, uJ)

        end
    end

    ### compute outgoing flow
    γout = self.Jmat * γin

    return [γin, γout]
end


function LP_1x1!(γin::Vector{T}, γin_max::Vector{R}, γout_max::Vector{S}, uJ::Vector{P}) where {T,R,S,P <: Real}
    γin[1] = min( γin_max[1], (1-uJ[1]) * γout_max[1] ) 
    return nothing
end

function LP_1x2!(self::Junction,
                γin::Vector{T}, γin_max::Vector{R}, γout_max::Vector{S}, uJ::Vector{P}; 
                tol=Float64x2(1e-10)) where {T,R,S,P <: Real}

    #TODO: que faire dans le cas u=0 et γmax=0 ?
    ### build controlled simplex Ω(u)
    α = self.α[1,1] #? self.α est constant
    γ1L_max = (1-uJ[1]+tol) * γout_max[1] 
    γ2L_max = (1-uJ[2]+tol) * γout_max[2]

    ### adjust repartition matrix to the control
    x = uJ[1] - uJ[2]
    alpha = 0.5*x*(x-1) + α*(1-x^2) + x*tol^2
    self.Jmat[1,1] = alpha #? self.Jmat dépend de u
    self.Jmat[2,1] = 1-alpha

    ### solve LP
    γin[1] = min( γin_max[1], min( γ1L_max/alpha, γ2L_max/(1-alpha) ) )

    return nothing
end

function LP_2x1!(self::Junction,
                γin::Vector{T}, γin_max::Vector{R}, γout_max::Vector{S}, 
                uJ::Vector{P}) where {T,R,S,P <: Real}

    ### build controlled simplex Ω(u)
    q = self.priorities[1]
    γ1R, γ2R = 0., 0.
    γ1R_max = γin_max[1]
    γ2R_max = γin_max[2]
    γ3L_max = (1 - uJ[1]) * γout_max[1] 
    
    ### solve LP
    if γ1R_max + γ2R_max <= γ3L_max
        γ1R = γ1R_max
        γ2R = γ2R_max
    else
        if γ1R_max >= q*γ3L_max && γ2R_max >= (1-q)*γ3L_max
            γ1R = q*γ3L_max
            γ2R = (1-q)*γ3L_max
        elseif γ1R_max < q*γ3L_max && γ2R_max >= (1-q)*γ3L_max
            γ1R = γ1R_max
            γ2R = γ3L_max - γ1R_max
        elseif γ1R_max >= q*γ3L_max && γ2R_max < (1-q)*γ3L_max
            γ1R = γ3L_max - γ2R_max
            γ2R = γ2R_max 
        end
    end

    γin[1] = γ1R 
    γin[2] = γ2R

    return nothing
end

function LP_2x2!(self::Junction,
    γin::Vector{T}, γin_max::Vector{R}, γout_max::Vector{S}, uJ::Vector{P}; 
    tol=Float64x2(1e-10)) where {T,R,S,P <: Real}

    tol *= 0.
    # println("---------- 2x2 ----------")
    ### build controlled simplex Ω(u)
    #? self.α est constant
    α = self.α[1,1]
    β = self.α[1,2] 
    γ1R_max = γin_max[1] 
    γ2R_max = γin_max[2]
    γ3L_max = (1 - uJ[1]) * γout_max[1];
    γ4L_max = (1 - uJ[2]) * γout_max[2]; 

    ### adjust repartition matrix to the control
    x = uJ[1] - uJ[2]
    alpha = 0.5*x*(x-1) + α*(1-x^2) + x*tol^2
    beta  = 0.5*x*(x-1) + β*(1-x^2) + x*tol^2
    self.Jmat[1,1] = alpha   ; self.Jmat[1,2] = beta 
    self.Jmat[2,1] = 1-alpha ; self.Jmat[2,2] = 1-beta
    # @show alpha, beta
    ### create buffer 
    γ1R, γ2R = 0., 0.


    ### check if "2x1-like" or genuine 2x2 
    if uJ[1] - uJ[2] == -1
        self.priorities = [0.5; 0.5]
        LP_2x1!(self, γin, γin_max, γout_max[1,:], uJ)
    elseif uJ[1] - uJ[2] == 1
        self.priorities = [0.5; 0.5]
        LP_2x1!(self, γin, γin_max, γout_max[2,:], uJ)
    else # genuine 2x2    
        Δ = alpha*(1-beta) - beta*(1-alpha)
        Γ_1 = (  (1-beta)*γ3L_max - beta*γ4L_max ) / Δ
        Γ_2 = ( -(1-alpha)*γ3L_max + alpha*γ4L_max ) / Δ
        if Γ_1 <= γ1R_max && Γ_2 <= γ2R_max
            # println("----- case 1 -----")
            if Γ_1 >= 0 && Γ_2 < 0
                γ2R = 0
                if alpha < beta 
                    γ1R = γ3L_max / alpha
                else
                    γ1R = γ4L_max / (1-alpha)
                end
            elseif Γ_1 < 0 && Γ_2 >= 0
                γ1R = 0
                if alpha < beta
                    γ2R = γ3L_max / beta
                else
                    γ2R = γ4L_max / (1-beta)
                end
            elseif Γ_1 < 0 && Γ_2 < 0 
                γ1R = 0
                γ2R = 0
            else
                γ1R = Γ_1
                γ2R = Γ_2
            end
        elseif Γ_1 > γ1R_max && Γ_2 > γ2R_max
            # println("----- case 2 -----")
            γ1R = γ1R_max
            γ2R = γ2R_max
        elseif Γ_1 > γ1R_max && Γ_2 <= γ2R_max
            # println("----- case 3 -----")
            if alpha < beta
                # println("----- case 3.1 -----")
                if γ3L_max-alpha*γ1R_max >= 0
                    # println("--- case 3.1.1 ---")
                    γ1R = γ1R_max
                    γ2R = min( (γ3L_max - alpha*γ1R_max) / beta, γ2R_max )
                else
                    # println("--- case 3.1.2 ---")
                    γ1R = γ3L_max/alpha
                    γ2R = 0.
                end
            else
                # println("----- case 3.2 -----")
                if γ4L_max-(1-alpha)*γ1R_max >= 0
                    # println("--- case 3.2.1 ---")
                    γ1R = γ1R_max
                    γ2R = min( (γ4L_max-(1-alpha)*γ1R_max) / (1-beta), γ2R_max)
                else
                    # println("--- case 3.2.2 ---")
                    γ1R = γ4L_max/(1-alpha)
                    γ2R = 0.
                end
                # 
            end
        elseif Γ_1 <= γ1R_max && Γ_2 > γ2R_max
            # println("----- case 4 -----")
            if alpha > beta
                # println("----- case 4.1 -----")
                if γ3L_max - beta*γ2R_max >= 0
                    # println("--- case 4.1.1 ---")
                    γ1R = min( (γ3L_max-beta*γ2R_max)/alpha, γ1R_max)
                    γ2R = γ2R_max
                else
                    # println("--- case 4.1.2 ---")
                    γ1R = 0
                    γ2R = γ3L_max/beta 
                end
            else
                # println("----- case 4.2 -----") 
                if γ4L_max - (1-beta)*γ2R_max >= 0
                    # println("--- case 4.2.1 ---")
                    γ1R = min( (γ4L_max - (1-beta)*γ2R_max)/(1-alpha), γ1R_max)
                    γ2R = γ2R_max
                else
                    # println("--- case 4.2.2 ---")
                    γ1R = 0
                    γ2R = γ4L_max / (1-beta)
                end
            end
        end
    end

    γin[1] = γ1R
    γin[2] = γ2R

    return nothing
end 