export Model, init, ComputePrimal!, AutoJac!

struct Mesh
    L               :: Float64 # length of one road 
    NumberOfCenters :: Int64   # mesh size of one road
    Δx              :: Float64 # area of one cell
    nodes           :: Vector{Float64}
    centers         :: Vector{Float64}

    function Mesh(L, Nc, Nr) 
        Δx = L/Nc
        nodes = [j*Δx for j in 0:Nr*Nc]
        centers = [(nodes[j]+nodes[j+1])/2 for j in 1:Nr*Nc]
        new(L, Nc, Δx, nodes, centers)
    end
end

# controlled macroscopic model of traffic flow based on LWR
struct Model
    net          :: Network                  # A road network in a weighted graph structure
    mapping      :: Vector{Vector{Int64}}    # mapping between Vectorial and Graph indices

    mesh         :: Mesh

    Tf           :: Float64                  # Final time
    Nt           :: Int64                    # Number of time step
    Δt_t         :: Vector{Float64}          # list of time steps (size = Nt - 1)
    Δt           :: Float64                  # Time step

    ρ_t          :: Vector{Vector{Float64}}        # densities over time
    γ_t          :: Vector{Vector{Float64}}        # Neumann bc at junctions over time
    u_t          :: Vector{Vector{Float64}}        # controls over time

    FV           :: Vector{Float64}                # FV flow (at a given time: temporary data)
    γ            :: Vector{Float64}                # Neumann bc (at a given time: temporary data)

    M_FV_t       :: Vector{Matrix{Float64}}       # table of FV jacobians over time
    M_LP_t       :: Vector{Matrix{Float64}}       # table of LP jacobians over time
    N_t          :: Vector{Matrix{Float64}}       # table of matrices N over time

    diffresult_FV :: Matrix{Float64}
    diffresult_LP :: Matrix{Float64}
end

#*** CONSTRUCTORS ***
function Model(net::Network, Nc; cfl=0.5, Tf=1.0)
    L=1.0  
    Nr = net.NumberOfRoads

    mesh = Mesh(L, Nc, Nr)

    Δt = cfl*mesh.Δx/maximum(net.vmax_t) 
    Nt = Int64(floor(Tf/Δt + 1)) 
    Δt_t = zeros(Nt-1)
    
    ρ_t = Vector{Vector{Float64}}(undef, Nt) 
    u_t = Vector{Vector{Float64}}(undef, Nt-1)
    γ_t = similar(u_t) 

    mapping = construct_mapping(net, mesh)

    FV = Vector{Float64}(undef, Nr*Nc) # zeros(Nr*Nc) 
    γ  = Vector{Float64}(undef, 2*Nr)  # zeros(2*Nr) 

    M_FV_t = Vector{Matrix{Float64}}(undef, Nt-1)
    M_LP_t = similar(M_FV_t) # [zeros(Nr*Nc, Nr*Nc) for _ in 1:Nt]
    N_t    = similar(M_FV_t) # [zeros(Nr*Nc, Nr) for _ in 1:Nt]

    diffresult_FV = Matrix{Float64}(undef, Nr*Nc, Nr*(Nc+2))
    diffresult_LP = Matrix{Float64}(undef, 2*Nr, Nr*(Nc+1))

    Model(net, mapping, mesh, Tf, Nt, Δt_t, Δt, ρ_t, γ_t, u_t,
    FV, γ, M_FV_t, M_LP_t, N_t, diffresult_FV, diffresult_LP)
end

#*** INITIALIZERS ***
function init(model::Model)
    Nr = model.net.NumberOfRoads; Nc = model.mesh.NumberOfCenters
    for n in 1:model.Nt-1
        model.ρ_t[n] = zeros(Nr*Nc)
        model.γ_t[n] = zeros(2*Nr)
        model.u_t[n] = zeros(Nr)

        model.M_FV_t[n] = zeros(Nr*Nc, Nr*Nc)
        model.M_LP_t[n] = zeros(Nr*Nc, Nr*Nc)
        model.N_t[n] = zeros(Nr*Nc, Nr)
    end
    model.ρ_t[end] = zeros(Nr*Nc)

    fill!(model.FV, 0)
    fill!(model.γ, 0)
    fill!(model.diffresult_FV, 0)
    fill!(model.diffresult_LP, 0)

    nothing
end

#*** GETTERS *** 
### indices ordering for each road : ρ1, ..., ρNc, γL, γR

function get_global_index_field(model::Model, i_edge, i_cell)
    return i_cell != :end ? model.mapping[i_edge][i_cell] : model.mapping[i_edge][end-2]
end

function get_global_index_flux(model::Model, i_edge, side)
# Nc + 1 --> γL
# Nc + 2 --> γR
    return side < 0 ? model.mapping[i_edge][end-1] : model.mapping[i_edge][end]
end

#*** SETTERS ***
function set_control!(model::Model, v_t::Vector{Vector{Float64}})
    model.u_t .= copy(v_t)
end
    

#*** METHODS ***
function construct_mapping(net::Network, mesh::Mesh)
    mapping = Vector{Vector{Int64}}(undef, net.NumberOfRoads)
    global_idx = 1
    @inbounds for roadindex = 1:net.NumberOfRoads
        mapping[roadindex] = zeros(Int64, mesh.NumberOfCenters+2)
        @inbounds for j = 1:mesh.NumberOfCenters
            mapping[roadindex][j] = global_idx
            global_idx += 1
        end
        mapping[roadindex][mesh.NumberOfCenters+1] = 2*roadindex - 1
        mapping[roadindex][mesh.NumberOfCenters+2] = 2*roadindex
    end
    return mapping 
end


####### Differentiables functions ####### 
function f_FV!(FV, model::Model, ρ, γ)  


    #! compatibilité ave ForwardDiff.Dual
    FV = get_tmp(FV, [ρ; γ])

    @inbounds for k in 1:model.net.NumberOfRoads # road k

        j = 1 #? left cell

        id_γL = get_global_index_flux(model, k, -1)
        id_ρC = get_global_index_field(model, k, j)
        id_ρR = get_global_index_field(model, k, j+1)
        Fp = Fnum(model.net, ρ[id_ρC], ρ[id_ρR], k) 
        Fm = γ[id_γL] 

        FV[id_ρC] = (-1/model.mesh.Δx)*(Fp - Fm)  

        @inbounds for j = 2:model.mesh.NumberOfCenters-1 #? central cells

            id_ρL = get_global_index_field(model, k, j-1)
            id_ρC = get_global_index_field(model, k, j)
            id_ρR = get_global_index_field(model, k, j+1)
            Fp = Fnum(model.net, ρ[id_ρC], ρ[id_ρR], k) 
            Fm = Fnum(model.net, ρ[id_ρL], ρ[id_ρC], k) 

            FV[id_ρC] = (-1/model.mesh.Δx)*(Fp - Fm)
            
        end

        j = model.mesh.NumberOfCenters #? right cell

        id_ρL = get_global_index_field(model, k, j-1)
        id_ρC = get_global_index_field(model, k, j)
        id_γR = get_global_index_flux(model, k, 1)
        Fp = γ[id_γR] 
        Fm = Fnum(model.net, ρ[id_ρL], ρ[id_ρC], k) 

        FV[id_ρC] = (-1/model.mesh.Δx)*(Fp - Fm) 

    end

    return FV
end


function ϕ!(γ, model::Model, ρ, u) 

    #! compatibilité ave ForwardDiff.Dual
    γ = get_tmp(γ, [ρ; u])
    
    g = model.net.graph

    @inbounds for v in vertices(g)

        Inc = inneighbors(g, v)
        Out = outneighbors(g, v)

        #TODOOO
        # Inc = model.net.Junctions[v].InVertices
        # Out = model.net.Junctions[v].OutVertices

        Inc_roads = [model.net.Adj_mat[idx_in, v] for idx_in in Inc]
        Out_roads = [model.net.Adj_mat[v, idx_out] for idx_out in Out]

        #! preallocate (dualcache ?)
        γin_max = @views model.net.Junctions[v].γin_max
        γout_max = @views model.net.Junctions[v].γout_max


        if length(Inc)*length(Out) > 0 #? the node is not a leaf of the graph

            #? Get max ingoing fluxes
            @inbounds for k in eachindex(Inc_roads) #(k, idx_in) in enumerate(Inc_roads)
                #? get the global (vectorial) coordinates of the field
                id_field = get_global_index_field(model, Inc_roads[k], :end)
                #? compute and store maximal flux
                γin_max[k] = fmaxR(model.net, ρ[id_field], Inc_roads[k])
            end

            #? Get max outgoing fluxes
            @inbounds for k in eachindex(Out_roads) #(k, idx_out) in enumerate(Out_roads)
                #? get the global (vectorial) coordinates of the field
                id_field = get_global_index_field(model, Out_roads[k], 1)
                #? compute and store maximal flux (this one is controlled)
                γout_max[k] = fmaxL(model.net, ρ[id_field], Out_roads[k])               
            end

            #? Solve LP problem
            γin, γout = LP(model.net.Junctions[v], γin_max, γout_max, u[Out_roads]) #uJ

            #? Update bc γR
            @inbounds for k in eachindex(Inc_roads) 
                #? get the global (vectorial) coordinates of the flow (bc)
                id_flux = get_global_index_flux(model, Inc_roads[k], 1)
                #? fill vectorial flux γ
                γ[id_flux] = γin[k]
            end

            #? Update bc γL
            @inbounds for k in eachindex(Out_roads) 
                #? get the global (vectorial) coordinates of the flow (bc)
                id_flux = get_global_index_flux(model, Out_roads[k], -1)
                #? fill vectorial flux γ
                γ[id_flux] = γout[k]
            end

        else #? the node is a leaf (thus this is a 1x1 junction)
            
            if length(Inc) == 0  #? leaf is an input of the graph (γL)

                #? Get max ingoing flux
                index = 0
                @inbounds for k in eachindex(model.net.list_leaves) #(k,l) in enumerate(model.net.list_leaves)
                    if v == model.net.list_leaves[k] #l
                        index = k
                        break 
                    end
                end

                γin_max[1]  = model.net.Ψ[index]

                #? Get max outgoing flux
                idx_out = Out_roads[1]
                id_field = get_global_index_field(model, idx_out, 1)
                
                #! no control here (reactivate for single junction debug)
                γout_max[1] = fmaxL(model.net, ρ[id_field], idx_out)
                
                #? Solve LP problem
                γin, γout = LP(model.net.Junctions[v], γin_max, γout_max, u[Out_roads]) # uJ
                
                #? Update bc γL
                id_flux = get_global_index_flux(model, idx_out, -1)
                γ[id_flux] = γout[1]
                
            else  #? leaf is an output of the graph (γR)
            
                #? Get max ingoing flux
                idx_in = Inc_roads[1]
                id_field = get_global_index_field(model, idx_in, :end)
                γin_max[1] = fmaxR(model.net, ρ[id_field], idx_in)
                
                #? Get max outgoing flux
                index = 0
                @inbounds for k in eachindex(model.net.list_leaves) #(k,l) in enumerate(model.net.list_leaves)
                    if v == model.net.list_leaves[k] #l
                        index = k
                        break 
                    end
                end

                #TODOO : mettre une grande valeur ? (pas de contrainte sur le flux sortant d'une feuille du graphe)
                #TODOO : comparer avec le fait de mettre γ[id_flux] = f(ρR)
                # γout_max[1] = γin_max[1] #no restriction #model.net.Ψ[index]
                γout_max[1] = 0.25 #model.net.Ψ[index] #! threeways_graph
                # γout_max[1] = model.net.Ψ[index]

                #? Solve LP problem
                γin, γout = LP(model.net.Junctions[v], γin_max, γout_max, [0.])

                #? Update bc γR
                id_flux = get_global_index_flux(model, idx_in, 1)
                γ[id_flux] = γin[1]
                
            end
        end
    end

    return γ
end


"""
    flot_FV(model::Model, ρ, γ)

construct the FV flow (ODE).
"""
flot_FV!(model::Model, ρ, γ) = f_FV!(model.FV, model, ρ, γ)

neumann_bc!(model::Model, ρ::Vector{Float64}, u::Vector{Float64}) = ϕ!(model.γ, model, ρ, u)


# function my_get_tmp(dc::PreallocationTools.DiffCache, u::Vector{T}) where {T <: ForwardDiff.Dual}
#     return reshape(reinterpret(T, dc.dual_du), size(dc.du))
# end
# function my_get_tmp(dc::PreallocationTools.DiffCache, u::Vector{T}) where {T <: Float64}
#     return dc.du
# end

"""
    AutoJac!(model::Model, n::Int64)

Jacobian computed with ForwardDiff without allocations
"""
function AutoJac!(model::Model, n::Int64)

    Nc = model.mesh.NumberOfCenters
    Nr = model.net.NumberOfRoads

    foo!(FV, rg) = f_FV!(FV, model, view(rg, 1:Nc*Nr), view(rg, Nc*Nr + 1:Nr*Nc + 2Nr))
    bar!(LP, ru) = ϕ!(LP, model, view(ru, 1:Nc*Nr), view(ru, Nc*Nr + 1 : Nr*Nc + Nr) )

    ### mask control ### 
    # model.u_t[n][1] = 0.

    rg_n = [model.ρ_t[n]; model.γ_t[n]]
    ru_n = [model.ρ_t[n]; model.u_t[n]]

    chunksize_rg = Val{ForwardDiff.pickchunksize(length(rg_n))}.parameters[1]
    chunksize_ru = Val{ForwardDiff.pickchunksize(length(ru_n))}.parameters[1]

    cfg_fv = ForwardDiff.JacobianConfig(foo!, model.FV, rg_n, ForwardDiff.Chunk{chunksize_rg}()) # model.net.Ne * (Nc+2)
    cfg_lp = ForwardDiff.JacobianConfig(bar!, model.γ, ru_n, ForwardDiff.Chunk{chunksize_ru}()) # model.net.Ne * (Nc+1)

    ForwardDiff.jacobian!(model.diffresult_FV, foo!, model.FV, rg_n, cfg_fv)
    ForwardDiff.jacobian!(model.diffresult_LP, bar!, model.γ, ru_n, cfg_lp)

    #! we assume that all the roads have the same number of cells
    dγFV = model.diffresult_FV[:, Nc*Nr+1 : end]
    drϕ  = model.diffresult_LP[:, 1:Nc*Nr]
    duϕ  = model.diffresult_LP[:, Nc*Nr+1 : end]

    model.M_FV_t[n] .= model.diffresult_FV[:, 1:Nc*Nr] # --> drFV
    mul!(model.M_LP_t[n], dγFV, drϕ)
    mul!(model.N_t[n], dγFV, duϕ)

    ## mask control ### 
    # [mat[:,1] .= 0. for mat in model.N_t]

    nothing
end


function ComputePrimal!(model::Model, initialstate; computejacobians=false)

    # initialization 
    model.ρ_t[1] .= initialstate(model.mesh.centers)
    neumann_bc!(model, model.ρ_t[1], model.u_t[1])
    model.γ_t[1] .=  model.γ
    computejacobians && AutoJac!(model, 1)

    time = 0.
    dt   = model.Δt

    @inbounds for n in 2:model.Nt
        if time + dt > model.Tf
            dt = model.Tf - time 
        end

        # transport step 
        flot_FV!(model, model.ρ_t[n-1], model.γ_t[n-1])
        RK1!(model.ρ_t[n], model.ρ_t[n-1], model.FV, dt)

        if n<model.Nt 
            ### boundary conditions update 
            neumann_bc!(model, model.ρ_t[n], model.u_t[n])
            model.γ_t[n] .= model.γ

            # autodiff
            computejacobians && AutoJac!(model, n) 
        end 
        # updates 
        time += dt 
        model.Δt_t[n-1] = dt
    end

    nothing 
end
