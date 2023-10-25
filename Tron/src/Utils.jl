export init_dir, wopen, RM, load_data
export Pol!, aν, mν, Mν, d1Mν, d2Mν
export d1mν, d2mν, ∇TVmat
export spot_Λ
export l2_norm
export extend_lane, extend_rho!, extend_u!
export constant_J, constant_gamma
export write_files, retrieve_data
export dict2plot
export get_global_index_field
export test_Optimizer

### deep-copy of struct
Base.copy(x::T) where T = T([getfield(x, k) for k ∈ fieldnames(T)]...)

### global variables ### #! vérifier si ça ne pose pas de problèmes 
global _NU = 1e-16
# global _EPS = 1e7
# global _ALPHA = 5
# global VERBOSE = false
# global tol = 1e-5

### system ###

# removes a directory with all its subdirectories
RM(x) = rm(x, force=true, recursive=true)

# creates and moves to an empty directory
function init_dir(dir)
    # creates ./dir/ directory if it doesn't exist
    isdir(dir) || mkdir(dir)
    # moves to the directory
    cd(dir); println("moving to ", pwd())
    # cleans the directory if not empty
    RM.(readdir())
end

# opens and writes an array into a file with specified opening mode
function wopen(filename, openmode, array)
    open(filename, openmode) do f
        writedlm(f, array)
    end
end

### maths ###

##### Lagrange interpolation at junctions 1x2 
function Lagrange(args)
    x,α = args 
    return 0.5*x*(x-1) + α*(1-x^2) + x*tol^2
end

function Pol!(args)
    Junc, u1, u2 = args  
    alpha = Lagrange([u1-u2, Junc.Jmat[1,1]])
    Junc.Jmat[1,1] = alpha 
    Junc.Jmat[2,1] = 1 - alpha 
end


##### LogSumExp (LSE) #####
function LSE(X; α=50., tol=1e-15)
    xmax = maximum(X)
    kmax = argmax(X)
    N = length(X)
    tmp = 0
    for i = 1:N
        if i != kmax
            tmp += exp(α*(X[i]- xmax + tol))
        end
    end
    # return log(tmp-(length(X)-1))/α
    return log(tmp)/α
end

minLSE(X; α=10.) = - LSE(-X; α=α)
maxLSE(X; α=10.) =   LSE(X ; α=α)

##### SOFTMAX #####
function Softmax(X; α=50.) #TODO : this version should avoid overflow
    tmp1, tmp2 = 0,0
    xmax = maximum(X)
    kmax = argmax(X)
    for i in eachindex(X)
        if i != kmax
            eax = exp(α*(X[i]-xmax))
            tmp1 += X[i]*eax
            tmp2 += eax
        end
    end
    return (xmax + tmp1) / (1 + tmp2)
end

function ∇Softmax(X; α=50.)
    n = length(X)
    vec = zeros(n)
    tmp = 0
    S = Softmax(X, α=α) 
    for j = 1:n
        eax = exp(α*X[j])
        tmp += eax
        vec[j] = eax*(1 + α*(X[j]-S))
    end
    return vec / tmp
end

minSM(x) = -Softmax(-x; α=  _ALPHA)
maxSM(x) =  Softmax( x; α=  _ALPHA) 

dminSM(x) = ∇Softmax(-x; α= _ALPHA)
dmaxSM(x) = ∇Softmax( x; α= _ALPHA) 

# regularized abs
function aν(v)
    return sqrt(v^2 + _NU^2)
end

# auxilliary function
function θν(args)
    x, y = args
    return (x-y) / aν(x-y)
end

# regularized minimum
function mν(args)
    x, y = args
    return 0.5*( x+y - aν(x-y)  )
    # return minimum(args) #! DEBUG
end
function d1mν(args)
    return 0.5*( 1 - θν(args) )
end
function d2mν(args)
    return 0.5*( 1 + θν(args) )
end

# regularized maximum
function Mν(args)
    x, y = args
    return 0.5*( x+y + aν(x-y)  )
    # return maximum(args) #! DEBUG
end
function d1Mν(args)
    return 0.5*( 1 + θν(args) )
end
function d2Mν(args)
    return 0.5*( 1 - θν(args) )
end

mν(x, y) = mν([x,y])
Mν(x, y) = Mν([x,y])
d1mν(x, y) = d1mν([x,y])
d2mν(x, y) = d2mν([x,y])
d1Mν(x, y) = d1Mν([x,y])
d2Mν(x, y) = d2Mν([x,y])
θν(x, y) = θν([x,y])

##### GENERIC EXTREMA #####  
regmin = min #minimum #minSM 
regmax = max #maximum #maxSM 


# projection on [a, b]^N
function proj(v; a=0., b=1.)
    return min.(max.(v, a), b)
end

function projν(v; a=0., b=1.)
    return mν.(Mν.(v, a), b)
end

# extract raws from a matrix
raws(M) = (view(M, i, :) for i in 1:size(M, 1))

# extract columns from a matrix
columns(M) = (view(M, :, j) for j in 1:size(M, 2))

# Total Variations
## of a vector
function TV(v)
    N = length(v)  
    tv = 0.
    for n = 1:N-1
        tv += aν(v[n+1]-v[n]) 
    end
    return tv
end
## of an array
function TVmat(M)
    tv = 0.
    for raw in raws(M)
        tv += TV(raw)
    end
    return tv #/size(M, 1) # normalization by number of roads
end
### of a list of vectors 
TotalVariations(list_vec) = TVmat(hcat(list_vec...))
    
function D(n,u)
    return u[n+1] - u[n]
end
function δ(i,j)
    return ifelse(i==j, 1, 0)
end

function ∇TVmat(M) 
    Nr = size(M, 1) # number of roads
    Nt = size(M, 2) # number of time steps
    A = zeros(Nr, Nt)
    sum_k = zeros(Nr, Nt)  
    sum_n = zeros(Nr, Nt)
    for (k,u) in enumerate(raws(M))
        for n = 1:Nt-1
            Dknu = D(n,u)
            A .*= 0.
            for m = 1:Nt-1
                for i = 1:Nr
                    A[i,m] = δ(k,i)*( δ(m,n+1) - δ(m,n) )
                end
            end
            sum_n .+= A * Dknu/aν(Dknu)
        end
        sum_k .+= sum_n
    end
    return sum_k
end

∇TotalVariations(list_vec) = columns( ∇TVmat(hcat(list_vec...)) )

function optimality_conditions(u, dJ; tol=0)
    for (k,uk) in enumerate(u)
        if uk == 0
            @test dJ[k] > -tol
        elseif uk == 1
            @test dJ[k] < tol
        elseif 0 < uk < 1
            @test -tol <= dJ[k] <= tol
        else
            @warn "encoutered invalid value in test : u must be between 0 and 1."
        end
    end
end 

function ComputeOptimalityConditions(u, dJ; M=1., λ=0.)
    Λ = copy(u) .* 0
    for n in eachindex(Λ) # loop over time 
        for r in eachindex( Λ[n] ) # loop over roads 
            Λ[n][r] = proj(dJ[n][r] + λ, a=u[n][r] - M, b=u[n][r]) 
        end
    end
    return norm(Λ)
end

# function spot_Λ(u,dJ,Λ)
#     SΛ = zeros(Int8, size(Λ))
#     for n in 1:size(Λ)[2]
#         for k in 1:size(Λ)[1]
#             if Λ[k,n] == u[k,n]
#                 SΛ[k,n] += 3 
#             end
#             if Λ[k,n] == u[k,n]-1
#                 SΛ[k,n] += 10 
#             end
#             if Λ[k,n] == dJ[k,n]
#                 SΛ[k,n] += 17 
#             end
#             if SΛ[k,n] == 0
#                 @warn "SΛ[$k,$n] == 0"
#             end
#         end
#     end
#     return SΛ
# end

# -----------------------------------------------------------------------------


function l2_norm(M, weights)
    norm_M = 0.
    for nw in eachindex(weights)
        for nr in eachindex(M[1])
            norm_M += M[nw][nr]^2 * weights[nw]
        end
    end
    return norm_M
end

function extend_lane(lane, Ne::Int64, Nc::Int64)
    ext_lane = zeros(Ne*Nc)
    for k in 1:Ne
        ext_lane[(k-1)*Nc + 1 : k*Nc] .= lane[k]
    end
    return ext_lane
end


function extend_rho!(new_rho, old_rho)
    new_Nc = length(new_rho)
    old_Nc = length(old_rho)
    K = div(new_Nc, old_Nc)
    for j = 1:old_Nc
        new_rho[(j-1)*K + 1 : j*K] .= old_rho[j]
    end
end

function extend_u!(new_u, old_u)
    new_Nt = size(new_u,2)
    old_Nt = size(old_u,2)
    K = Int( ceil(new_Nt / old_Nt) )
    for j = 1:old_Nt
        new_u[:, (j-1)*K + 1 : min(new_Nt, j*K)] .= old_u[:, j]
    end
end

function constant_J(v, opt::Optim; idx=[1,2])                                                                                                                                                                                                  
    u= 0. * copy(opt.model.u_t)                                                                                                                                                                                             
    for el in u   
       el[idx] .= v                                                                                                                                                                                                            
    end                                                                                                                                                                                                                     
    ComputeCost!(opt, u=u)                                                                                                                                                                                                
end  


function constant_gamma(v, opt::Optim; idx=[3,4])                                                                                                                                                                                                  
    u= 0. * copy(opt.model.u_t)                                                                                                                                                                                             
    for el in u   
       el[idx] .= v                                                                                                                                                                                                            
    end                                                                                                                                                                                                                     
    ComputeCost!(opt, u=u)  
    return [opt.model.γ[4]; opt.model.γ[5]  ]                                                                                                                                                                                     
end  



# function that takes a model and writes phi, v, psi, grad_psi in files 
# (one directory per variable, one file per time step)
function write_files(opt::Optim; init_name="data", test_case="generic")
    
    init_dir(init_name)

    init_dir("PARAMETERS")

    open("test_case", "w") do file
        print(file, test_case)
    end 
    open("algo", "w") do file
        print(file, opt.algo)
    end 
    open("param.txt", "w") do file  
        println(file, "==== run infos ====")
        # println(file, "max iterations: ", opt.itermax) #todo
        # println(file, "initial control: ", opt.u_init) #todo
        # println(file, "convergence tol: ", opt.CV_tol) #todo
        # println(file, "max iter before forced GD: ", opt.maxiter_b4_gd) #todo
        println(file, " ")
        println(file, "==== optim infos ====")
        println(file, "final time: ", opt.model.Tf)
        println(file, "algorithm: ", opt.algo)
        println(file, "Nmax: ", opt.Nmax)
        println(file, "Staffing regularization (initial): ", opt.regularizations[1])
        println(file, "BV regularization (initial): ", opt.regularizations[2])
        println(file, "fp tol: ", opt.fp_tol)
        println(file, "fp trigger: ", opt.trigger)
        println(file, "fp trigger delay factor: ", opt.delay_factor)
        println(file, "total iterations: ", opt.iter)
        println(file, "residual at last iteration: ", opt.res)
        println(file, " ")
        println(file, "==== model and mesh infos ====")
        println(file, "length of each road: ", opt.model.mesh.L)
        println(file, "number of primal cells (centers) per road: ", opt.model.mesh.NumberOfCenters)
        println(file, "number of time steps: ", opt.model.Nt)
        println(file, " ")
        println(file, "==== network infos ====")
        println(file, "number of roads: ", opt.model.net.NumberOfRoads)
        println(file, "number of junctions: ", opt.model.net.NumberOfVertices)
        println(file, "number of junctions that are leaves of the graph: ", opt.model.net.NumberOfLeaves)
        println(file, "location of leaves: ", opt.model.net.list_leaves)
        println(file, "density values at leaves of the graph: ", opt.model.net.ρ_leaves_t)
        println(file, "maximal velocities at roads: ", opt.model.net.vmax_t)
        println(file, "maximal densities at roads: ", opt.model.net.ρmax_t)

    end

    cd("..") # back from: init_dir("PARAMETERS")

    wopen("redlane", "w", opt.redlane)
    wopen("adjmat", "w", opt.model.mapping)
    wopen("density_max", "w", opt.model.net.ρmax_t)
    wopen("NumberOfCenters", "w", opt.model.mesh.NumberOfCenters)
    wopen("NumberOfRoads", "w", opt.model.net.NumberOfRoads)
    wopen("NumberOfTimeSteps", "w", opt.model.Nt)
    wopen("final_time", "w", opt.model.Tf)
    ### save mesh 
    wopen("centers", "w", opt.model.mesh.centers)
    wopen("nodes", "w", opt.model.mesh.nodes)

    ### SAVES DIAGNOSTICS
    init_dir("DIAGNOSTICS")

    begin 
        ### save snapshots times 
        wopen("snaptimes", "w", opt.saving_times)
        wopen("snapindices", "w", opt.saving_indices)

        init_dir("bests")
        begin
            wopen("density_best", "w", opt.diagnostics.ρ_best)
            wopen("control_best", "w", opt.diagnostics.u_best)
        end
        cd("..") # back from: init_dir("bests")

        init_dir("lasts")
        begin
            wopen("density_last", "w", opt.diagnostics.ρ_last)
            wopen("control_last", "w", opt.diagnostics.u_last)
        end
        cd("..") # back from: init_dir("lasts")

        init_dir("histories")
        begin
            wopen("loss_history", "w", opt.diagnostics.losshistory)
            wopen("cost_history", "w", opt.diagnostics.costhistory)
            wopen("staff_history", "w", opt.diagnostics.staffhistory)
            wopen("totvar_history", "w", opt.diagnostics.totvarhistory)
            wopen("residual_history", "w", opt.diagnostics.residualhistory)
            wopen("optimality_history", "w", opt.diagnostics.optimalityhistory)
            if length(opt.diagnostics.descentstephistory) > 0
                wopen("descentstep_history", "w", opt.diagnostics.descentstephistory)
            end
            wopen("thetaStaff_history", "w", opt.diagnostics.thetaStaff_history)
            wopen("thetaBV_history", "w", opt.diagnostics.thetaBV_history)
        end
        cd("..") # back from: init_dir("histories")

    end

    cd("..") # back from: init_dir("DIAGNOSTICS")

    cd("..") # back from: init_dir(init_name)

    nothing
end


function retrieve_data(dir::String)
    D = Dict()
    for (root, dirs, file) ∈ walkdir(dir)
        if occursin("PARAMETERS", root)
            [D[fi] = read(joinpath(root, fi), String) for fi in file];
        else 
            [D[fi] = readdlm(joinpath(root, fi)) for fi in file];
        end
    end
    return D
end

function test_Optimizer(;final_time=0.3, Nc=5, Nmax=2, coef_pen=0, coef_reg=0, 
    itermax=3, CV_tol=1e-2, algo="GD",
    lane=Bool.([0]), initstate=Constant(value=0.66), fp_tol=0.0, fp_trigger=10, delay_factor=1.0,
    ismaskcontrol=false, u_init=0.0,
    network = graph_1x1(1.0), lr_init=1e-1, decay=1e-1
    )

    model = Model(network, Nc, Tf=final_time)
    init(model)

    optim = Optim(model, Nmax, coef_pen, coef_reg)
    set_scheduler(optim, lr_init=lr_init, decay=decay)
    set_redlane!(optim, lane)
    optim.fp_tol = fp_tol
    optim.trigger = fp_trigger
    optim.delay_factor = delay_factor
    optim.ismaskcontrol = ismaskcontrol
    [u .= u_init for u in optim.ukm]
    [u .= u_init for u in optim.uk]
    [u .= u_init for u in optim.ukp]

    optim(itermax=itermax, CV_tol=CV_tol, algo=algo, initstate=initstate)

    return optim

end

function control_embedding!(new_u_t, old_u_t)
    Nt_old = length(old_u_t)
    Nt_new = length(new_u_t)
    ratio = div(Nt_new, Nt_old)
    for id_time in 1:Nt_old 
        id_time
        id_begin = (id_time-1)*ratio + 1
        id_end = min(Nt_new, id_time*ratio)
        for id_time_loc in id_begin:id_end
            new_u_t[id_time_loc] .= copy(old_u_t[id_time])
        end
    end
end 

function gridsearch(;test_case="generic", final_time=0.3, Nc=5, Nmax=5, 
    itermax=3, CV_tol=1e-5, algo_list=["GD", "GDFP", "FP"],
    lane=Bool.([0]), 
    initstate=Constant(value=0.66), delay_factor=2,
    ismaskcontrol=false, u_init_list=[1],
    coef_pen_list=[1e-8], 
    coef_reg_list=[1e-8], 
    fp_tol_list=[1e-3],
    fp_trigger_list=[10]
    )

    test_Optimizer(); #trigger compilation
    init_dir("GRIDSEARCH")
    for algo in algo_list
        init_dir(algo)
        for u_init in u_init_list
            init_dir("u_init_"*string(u_init))
            for coef_pen in coef_pen_list
                coef_reg = coef_pen 
                for fp_tol in fp_tol_list
                    for fp_trigger in fp_trigger_list
                        println("-------------- NEW TESTCASE -------------")
                        @show u_init
                        @show algo 
                        @show coef_pen
                        @show coef_reg
                        @show fp_tol 
                        @show fp_trigger

                        optimize = test_Optimizer(final_time=final_time, Nc=Nc, Nmax=Nmax, 
                            coef_pen=coef_pen, coef_reg=coef_reg, itermax=itermax, CV_tol=CV_tol, 
                            algo=algo, lane=lane, initstate=initstate, fp_tol=fp_tol, 
                            fp_trigger=fp_trigger, delay_factor=delay_factor,
                            ismaskcontrol=ismaskcontrol, u_init=u_init)

                        diagnostics(optimize, test_case=test_case, 
                            dirname="diagnostics_"*algo*"_u_init_"*
                            string(u_init)*"_coef_pen_"*string(coef_pen)*
                            "_coef_reg_"*string(coef_reg)*"_fp_tol_"*
                            string(fp_tol)*"_fp_trigger_"*string(fp_trigger))
                        
                        optimize = nothing # does it reset memory ??
                        GC.gc(true)
                        # if algo != "GDFP"
                        #     break # skip loop over fp_trigger except for GDFP
                        # end
                    end 
                end 
            end 
            cd("..")
        end
        cd("..")
    end 
    cd("..")
end


function dict2plot(D; dirname="dict2plot", sizegraph=(1200,800),
    loss=true, cost=true, staff=true, totvar=true)

    init_dir(dirname)

    init_dir("snapshots")
    begin
        for (kt, t_snap) ∈ enumerate(D["snaptimes"])
            println("kt = $kt, t_snap = $t_snap")
            dict2plot_density(D, kt, t_snap, iter2plot="best") 
            dict2plot_density(D, kt, t_snap, iter2plot="last")
        end   
    end
    cd("..") # exit from init_dir("snapshots")
    println("moving to ", pwd())

    init_dir("controls")
    begin 
        dict2plot_control(D, iter2plot="best")
        dict2plot_control(D, iter2plot="last")
    end
    cd("..") # exit from init_dir("controls")
    println("moving to ", pwd())

    init_dir("histories")
    begin 
        #TOODO: diagnostics.foo_history, diagnostics.weighted_foo_history
        println("plot mean density...")
        dict2plot_meandensity(D, iter2plot="best") 
        dict2plot_meandensity(D, iter2plot="last")
        println(" done.")

        println("plot loss history...")
        dict2plot_losshistory(D, loss=loss, cost=cost, staff=staff, totvar=totvar) 
        println(" done.")

        println("plot residual history...")
        dict2plot_residualhistory(D) 
        println(" done.")

        println("plot optimality history...")
        dict2plot_optimalityhistory(D) 
        println(" done.")

        if D["algo"] != "FP"
            println("plot descent step history...")
            dict2plot_descentstephistory(D) 
            println(" done.")
        end
    end
    cd("..") # exit from init_dir("histories")
    println("moving to ", pwd())

    cd("..") # exit from init_dir(dirname)
    println("moving to ", pwd())
end

#TODOO: struct GraphTopology
function make_geometry_threeways()
    ### vertices are plotted with coordinates (xv[i], yv[i])
    xv = [0 5 10 15 20 25 30 35 40 45 13 18 28 38 8 13 23 33]
    yv = [0 0 0 0 0 0 0 0 0 0 3 3 3 3 -3 -3 -3 -3]
    ### edge number rd_num_t[i] is plotted with coordinates (A_t[i], B_t[i])
    A_t = []; B_t = []; rd_num_t = []
    # main road
    push!(A_t, [0 0], [5 0], [10 0], [15 0], [20 0], [25 0], [30 0], [35 0], [40 0])
    push!(B_t, [5 0], [10 0], [15 0], [20 0], [25 0], [30 0], [35 0], [40 0], [45 0])
    push!(rd_num_t, 1,2,4,6,7,8,9,10,11)
    # upper side
    push!(A_t, [10 0], [13 3], [18 3], [28 3], [38 3], [18 3], [28 3])
    push!(B_t, [13 3], [18 3], [28 3], [38 3], [40 0], [20 0], [30 0])
    push!(rd_num_t, 5,18,20,22,23,19,21)
    # lower side
    push!(A_t, [5 0], [8 -3], [13 -3], [23 -3], [13 -3], [23 -3], [33 -3])
    push!(B_t, [8 -3], [13 -3], [23 -3], [33 -3], [15 0], [25 0], [35 0])
    push!(rd_num_t, 3,12,14,16,13,14,17);

    return xv, yv, A_t, B_t, rd_num_t
end

function make_geometry_circle()
    ### vertices are plotted with coordinates (xv[i], yv[i])
    xv = [0 5 10 10 10 10 15 20]
    yv = [10 10 0 5 15 20 10 10]
    ### edge number rd_num_t[i] is plotted with coordinates (A_t[i], B_t[i])
    A_t = []; B_t = []; rd_num_t = []
    V1 = [9.99 0]; V2 = [10 5]; V3 = [15 10]; V4 = [10 15]; 
    V5 = [10.01 20]; V6 = [20 10]; V7 = [5 10]; V8 = [0 10];
    push!(A_t, V1, V2, V3, V3, V4, V5, V7, V7)
    push!(B_t, V2, V3, V4, V6, V7, V4, V2, V8)
    push!(rd_num_t, 1,2,3,4,5,6,7,8);

    return xv, yv, A_t, B_t, rd_num_t
end

function make_geometry_2x2()
    ### vertices are plotted with coordinates (xv[i], yv[i])
    xv = [0 0 5 10 10] 
    yv = [5 -5 0 5 -5]
    ### edge number rd_num_t[i] is plotted with coordinates (A_t[i], B_t[i])
    A_t = []; B_t = []; rd_num_t = []
    push!(A_t, [0 5], [0 -5], [5,0], [5,0]) 
    push!(B_t, [5 0], [5 0], [10,5], [10,-5]) 
    push!(rd_num_t, 1, 2, 3, 4);

    return xv, yv, A_t, B_t, rd_num_t
end

function make_geometry_2x1()
    ### vertices are plotted with coordinates (xv[i], yv[i])
    xv = [0 0 5 10]
    yv = [5 -5 0 0] 
    ### edge number rd_num_t[i] is plotted with coordinates (A_t[i], B_t[i])
    A_t = []; B_t = []; rd_num_t = []
    push!(A_t, [0 5], [0 -5], [5,0]) 
    push!(B_t, [5 0], [5 0], [10,0])
    push!(rd_num_t, 1, 2, 3);

    return xv, yv, A_t, B_t, rd_num_t
end

function make_geometry_1x2()
    ### vertices are plotted with coordinates (xv[i], yv[i])
    xv = [0 5 10 10] 
    yv = [0 0 5 -5]
    ### edge number rd_num_t[i] is plotted with coordinates (A_t[i], B_t[i])
    A_t = []; B_t = []; rd_num_t = []
    push!(A_t, [0 0], [5 0], [5 0])
    push!(B_t, [5 0], [10 5], [10 -5])
    push!(rd_num_t, 1, 2, 3);

    return xv, yv, A_t, B_t, rd_num_t
end

function make_geometry_1x1()
    ### vertices are plotted with coordinates (xv[i], yv[i])
    xv = [0 5 10]
    yv = [0 0 0]
    ### edge number rd_num_t[i] is plotted with coordinates (A_t[i], B_t[i])
    A_t = []; B_t = []; rd_num_t = []
    push!(A_t, [0 0], [5 0])
    push!(B_t, [5 0], [10 0])
    push!(rd_num_t, 1, 2);

    return xv, yv, A_t, B_t, rd_num_t
end

GraphGeometry = Dict()
GraphGeometry["threeways"] = make_geometry_threeways()
GraphGeometry["circle"] = make_geometry_circle()
GraphGeometry["2x2"] = make_geometry_2x2()
GraphGeometry["2x1"] = make_geometry_2x1()
GraphGeometry["1x2"] = make_geometry_1x2()
GraphGeometry["1x1"] = make_geometry_1x1()

function get_global_index_field(AdjMatrix::Array, i_edge, i_cell)
    return i_cell != :end ? Int(AdjMatrix[i_edge, i_cell]) : Int(AdjMatrix[i_edge, end-2])
end

println("load Utils.jl")