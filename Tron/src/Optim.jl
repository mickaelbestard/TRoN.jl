export Optim, init_Optim, set_redlane!, set_ctrl_lane
export ComputeDual!, ComputeCost! 
export diagnostics


Base.@kwdef mutable struct Scheduler  
    lr_init :: Float64
    decay :: Float64
end 

# exponential scheduler
function (s::Scheduler)(iteration)
    return s.lr_init*exp(-s.decay*iteration) #s.lr_init / (1 + s.decay*iteration)
end

# time-based scheduler
function (s::Scheduler)(iteration, lr, itermax)
    decay = s.lr_init/itermax
    return lr / (1 + decay*iteration)
end

#TODO: utiliser un dictionnaire pour stocker les paramètres pour que la struct puisse 
#TODO: être non mutable --> opt.optimdata
### Create and solve an optimal control problem
Base.@kwdef mutable struct Optim 
    ### problem parameters
    model           :: Model
    Nmax            :: Int64                    # Number of available controls
    p_t             :: Vector{Vector{Float64}}  # list of adjoint states

    redlane         :: Vector{Float64}          # mask of cells to consider in the functionnal
    ismaskcontrol   :: Bool                     # whether to control cells in the redlane 
    regularizations :: Vector{Float64}          # θ = (θ_S, θ_B)
    δ_t             :: Vector{Float64}          # step size over iterations
    trigger         :: Int64                    # interval between FP triggers in GDFP
    delay_factor    :: Float64                  # when true, multiply delays the next FP triggering at each use

    ### preallocated buffers for optimization iterations 
    uk              :: Vector{Vector{Float64}}
    ukp             :: Vector{Vector{Float64}}
    ukm             :: Vector{Vector{Float64}}

    J               :: Float64
    CT              :: Float64 
    Staffing        :: Float64 
    TotVar          :: Float64

    J_best          :: Float64
    CT_best         :: Float64 
    Staffing_best   :: Float64 
    TotVar_best     :: Float64

    dJ              :: Vector{Vector{Float64}}
    dCT             :: Vector{Vector{Float64}}
    dStaffing       :: Vector{Vector{Float64}}
    dTotVar         :: Vector{Vector{Float64}}

    res             :: Float64 # residual of control iterative method
    iter            :: Int64   # store the number of iterations during the optimization process
    diagnostics     :: Diagnostics # store all diagnostics
    saving_times    :: Vector{Float64} # times at which we want to plot
    saving_indices  :: Vector{Int64}

    algo            :: String

    fp_tol          :: Float64 # tolerance in FP algo

    scheduler :: Scheduler # learning rate scheduler

    function Optim(model::Model, Nmax=2, coef_pen=1e2, coef_reg=1e2, 
                    freq_savingtimes=0.1)
        
        Nr = model.net.NumberOfRoads
        Nc = model.mesh.NumberOfCenters
        Nt = model.Nt 

        p_t = 0. .* copy(model.ρ_t)
        
        J = CT = Staffing = TotVar = 0.
        J_best = CT_best = Staffing_best = TotVar_best = 0.
        
        dJ = [zeros(Nr) for _ in 1:Nt-1]
        dCT = [zeros(Nr) for _ in 1:Nt-1]
        dStaffing = [zeros(Nr) for _ in 1:Nt-1]
        dTotVar = [zeros(Nr) for _ in 1:Nt-1]

        redlane = zeros(Bool, Nr*Nc)
        ismaskcontrol = false

        regularizations = [coef_pen, coef_reg]

        trigger = 5
        delay_factor = 1.0
        δ_t = []

        uk = 0. .* copy(model.u_t)
        ukp = 0. .* copy(model.u_t)
        ukm = 0. .* copy(model.u_t) 

        res = 0.
        iter = 0

        diagnostics = Diagnostics(model)

        nb_savingtimes = Int( floor( freq_savingtimes * Nt ) )
        nb_savingtimes = max(2, nb_savingtimes)
        saving_times = collect( range(0.0, stop=model.Tf, length= nb_savingtimes) )
        saving_indices = Int.( floor.( collect(range(1, Nt, nb_savingtimes)) ) )

        algo = "GD"
        fp_tol = 0.

        scheduler = Scheduler(lr_init=1e-1, decay=1e-1)

        new(model, Nmax, p_t, redlane, ismaskcontrol, regularizations, δ_t, 
        trigger, delay_factor, uk, ukp, ukm, J, CT, Staffing, TotVar,
        J_best, CT_best, Staffing_best, TotVar_best, 
        dJ, dCT, dStaffing, dTotVar, 
        res, iter, diagnostics, 
        saving_times, saving_indices, algo, fp_tol, scheduler)
    end
end # class Optim

function init_Optim(opt::Optim) 
    #TODO
end

#*** SETTERS ***

function set_scheduler(opt::Optim; lr_init=1e-1, decay=1e-1)
    @show opt.scheduler.lr_init = lr_init  
    @show opt.scheduler.decay = decay
end

function set_penalization(opt::Optim, value) 
    opt.regularizations[1] = value
end

function set_regularization(opt::Optim, value) 
    opt.regularizations[2] = value
end

"""
    set_regularizations(opt, index, value)

fill vector of regularization parameters θ = (θs, θb) at index `index`
"""
function set_regularizations(opt, index, value)
    opt.regularizations[index] = value
end

function set_ismaskcontrol(opt::Optim, boolean) 
    opt.ismaskcontrol = boolean
end

function set_redlane!(opt::Optim, lane)
    opt.redlane .= 0.
    @inbounds for k in 1:opt.model.net.NumberOfRoads
        if lane[k] #check if k-th road is in redlane
            id_ρL = get_global_index_field(opt.model, k, 1)
            id_ρR = get_global_index_field(opt.model, k, :end)
            opt.redlane[id_ρL : id_ρR] .= 1
        end
    end
end

#### NB : dans certains cas-tests (jonctions isolées 1 x m) on ne contrôle pas la 
### feuille entrante, d'où :
### mask = 1 .- redlane, suivi de
### mask[1] = 1 
function maskcontrol!(opt::Optim, X=opt.u_t, mask= 1 .- opt.redlane) 
    if ismaskcontrol
        println("=== ismaskcontrol ===")
        Nr = opt.model.net.NumberOfRoads
        Nc = opt.model.mesh.NumberOfCenters
        for x in X
            roadindex = 1
            @inbounds for j = 1:Nc:Nr*Nc 
                x[roadindex] .*= mask[j]
                roadindex += 1
            end 
        end
    end
end

function no_control_lane!(opt::Optim, X=opt.model.u_t, mask=opt.redlane) 
    if opt.ismaskcontrol
        # println("@@@@@@@@@@@@ no_control_lane!()")
        for x in X
            x[1] *= 0. #? on ne contrôle pas la première route (sauf cas 1xm)
        end
        # Nr = opt.model.net.NumberOfRoads
        # Nc = opt.model.mesh.NumberOfCenters
        # k=1
        # @inbounds for j = 1:Nc:length(opt.redlane) 
        #     X[k] .*= 1-opt.redlane[j]
        #     k += 1
        # end 
    end
    return X
end


function ComputeDual!(opt::Optim, initdual=opt.redlane) 
   
    ### reverse time (now 0 -> t -> T)
    reverse!(opt.p_t)#, dims=2) 
    reverse!(opt.model.Δt_t) 
    reverse!(opt.model.M_FV_t)
    reverse!(opt.model.M_LP_t)

    ### initial (terminal) condition
    opt.p_t[1] = initdual

    # transpose all matrices before loop
    Z = transpose(opt.model.M_FV_t + opt.model.M_LP_t)
    
    @inbounds for n in 2:opt.model.Nt
        opt.p_t[n] .= (I + opt.model.Δt_t[n-1] .* Z[n-1]) * opt.p_t[n-1] 
    end
   
    ### reverse time back (now T -> t -> 0)
    reverse!(opt.model.Δt_t) 
    reverse!(opt.model.M_FV_t)
    reverse!(opt.model.M_LP_t)
    reverse!(opt.p_t) 

    nothing 
end

#! faire attention à quelle instance d'Optim on donne en param!
function ComputeCost!(opt::Optim; u=opt.model.u_t, initialfunc = Constant(value=0.66)) 
    
    θs, θb = opt.regularizations

    set_control!(opt, u)
    ComputePrimal!(opt.model, initialfunc) 

    opt.CT = dot(opt.redlane, opt.model.ρ_t[end])

    opt.Staffing *= 0.
    opt.TotVar *= 0.

    for n in eachindex(u)
        opt.Staffing += max(0., sum(u[n]) - opt.Nmax)^2 * opt.model.Δt_t[n]
    end

    opt.TotVar = TotalVariations(u)


    opt.J = opt.CT + θs/2 * opt.Staffing + θb * opt.TotVar 
end

function ComputeGradient!(opt::Optim, u=opt.model.u_t)
    θs, θb = opt.regularizations

    opt.dTotVar .= ∇TotalVariations(u)
    for n in eachindex(opt.dJ)
        opt.dCT[n] = transpose(opt.model.N_t[n]) * opt.p_t[n]
        opt.dStaffing[n] .= max(0., sum(u[n]) - opt.Nmax) #∀n, dSu[n] ∈ R^Nr
    end
    opt.dJ = opt.dCT + θs*opt.dStaffing + θb*opt.dTotVar
    ###TODO: "no control lane" is ambiguous: dJ_i=0 and/or u_i=0 ?
    no_control_lane!(opt, opt.dJ)

    return opt.dJ
end

set_control!(opt::Optim, v_t::Vector{Vector{Float64}}) = set_control!(opt.model, v_t)

function grad(opt::Optim; itermax=100)
    ### needed in 2x1 (and 1x1?) ###
    # if min(opt.J, opt.res) < 1e-1 
    #     opt.scheduler.lr_init = 1e-3
    # end 
    ###

    # @show descent_step = opt.scheduler(opt.iter)
    # lr = opt.iter==1 ? opt.scheduler.lr_init : opt.diagnostics.descentstephistory[end]
    # @show descent_step = opt.scheduler( opt.iter )
    @show descent_step = opt.scheduler.lr_init / 
                        (1 + opt.scheduler.decay*opt.iter)
    push!(opt.diagnostics.descentstephistory, descent_step)
    return descent_step * opt.dJ
end

function myprod(u,v)
    return [u[i] .* v[i] for i in eachindex(u)]
end

function fixed(opt::Optim; itermax=100)
    return opt.uk .- 
        (1. * map.(x -> x < -opt.fp_tol, opt.dJ) 
        + myprod(opt.uk, 1. * map.(x -> -opt.fp_tol <= x <= opt.fp_tol, opt.dJ) )
        )
end

function fixed_explore(opt::Optim)
    return opt.uk .- 1. * map.(x -> x < -opt.fp_tol, opt.dJ)
end

function hybrid(opt::Optim; itermax=100)
    @show opt.iter, opt.trigger, opt.iter % opt.trigger
    if opt.iter % opt.trigger == 0 #&& opt.CT > 1.0
        opt.trigger = Int(floor(opt.trigger*opt.delay_factor))
        println("next fp trigger at iteration: $(opt.trigger)")
        return fixed(opt, itermax=itermax)
    end
    return grad(opt, itermax=itermax)
end

function hybrid_explore(opt::Optim)
    # @show opt.iter, opt.trigger, opt.iter % opt.trigger
    if opt.iter % opt.trigger == 0 #&& opt.J > 1e-1
        opt.trigger = Int(floor(opt.trigger*opt.delay_factor))
        println("next fp_exp trigger at iteration: $(opt.trigger)")
        return fixed_explore(opt)
    end
    return grad(opt)
end

global first_time = true

function fixedthengrad(opt::Optim)
    if opt.iter > 200# opt.CT < 10
        println("GD iteration")
        if first_time
            opt.uk = opt.diagnostics.u_best
            opt.ukp = opt.diagnostics.u_best
            opt.model.u_t = u_best
            opt.model.ρ_t = opt.diagnostics.ρ_best
            global first_time=false
        end
        return grad(opt)
    end
    println("FP iteration")
    return fixed(opt)
end

DictOfOptimizers = Dict(
    "GD" => grad,
    "FP" => fixed,
    "GDFP" => hybrid,
    )


function (opt::Optim)(;itermax=3, CV_tol=1e-2, algo="GD", initstate=Constant(value=0.66))
    # opt.scheduler.decay = opt.scheduler.lr_init / itermax 
    opt.iter = 0
    opt.algo = algo
    optimizer = DictOfOptimizers[algo]
    opt.res = CV_tol + 1 # to enter while loop, won't be saved

    no_control_lane!(opt, opt.ukm)
    no_control_lane!(opt, opt.uk)
    no_control_lane!(opt, opt.ukp)

    ComputeCost!(opt, u=opt.uk)

    opt.J_best = opt.J 
    println("J$(opt.iter) = $(opt.J)")

    update_and_storage!(opt)

    #######
    opt.regularizations .= 1e-8
    #######
    while opt.iter <= itermax && opt.res > CV_tol
        println(">>>>> iteration $(opt.iter)") 

        no_control_lane!(opt, opt.ukm)
        no_control_lane!(opt, opt.uk)
        no_control_lane!(opt, opt.ukp)

        set_control!(opt, opt.uk)

        ComputePrimal!(opt.model, initstate, computejacobians=true)
        ComputeDual!(opt)
        ComputeGradient!(opt)

        descent = optimizer(opt, itermax=itermax)
        opt.ukp .= proj.( opt.uk .- descent )

        no_control_lane!(opt, opt.ukm)
        no_control_lane!(opt, opt.uk)
        no_control_lane!(opt, opt.ukp)

        ComputeCost!(opt, u=opt.ukp)
        
        ### Updates 
        opt.res = ComputeOptimalityConditions(opt.uk, opt.dJ)
        update_and_storage!(opt)

        println("J$(opt.iter-1) = $(opt.J), res$(opt.iter-1) = $(opt.res)")
        @show opt.CT, opt.Staffing, opt.TotVar
        @show norm(opt.uk .- opt.ukm)

        ########### CALLBACKS #############
        if opt.algo == "GDFP" && opt.res < 1 
            @show opt.regularizations[1] = 1e-4
            @show opt.regularizations[2] = 1e-6
        end

        
        if opt.algo == "FPGD" && opt.iter > 50 #200 #50 #200
            @show optimizer = DictOfOptimizers["GD"]
        end
        ##################################
    end
    ### End of optimization, store last computed values ###
    store_last!(opt)

    return opt 
end


"""
    update_and_storage!(opt::Optim)

TBW
"""
function update_and_storage!(opt::Optim; with_diagnostics=true)
    
    # update controls 
    opt.ukm .= copy(opt.uk)  
    opt.uk  .= copy(opt.ukp)

    if with_diagnostics
        push!(opt.diagnostics.losshistory, opt.J)
        push!(opt.diagnostics.costhistory, opt.CT)
        push!(opt.diagnostics.staffhistory, opt.Staffing)
        push!(opt.diagnostics.totvarhistory, opt.TotVar)

        if opt.iter > 0
            #! res = Λ, but this can be modified
            push!(opt.diagnostics.residualhistory, opt.res)
            push!(opt.diagnostics.optimalityhistory, opt.res) 
        end

        push!(opt.diagnostics.thetaStaff_history, opt.regularizations[1])
        push!(opt.diagnostics.thetaBV_history, opt.regularizations[2])

        # Store the value of ρ_t and u_t, when current iteration is optimal
        if opt.J < opt.J_best #&& opt.iter > 0 #! DEBUG
            println("### update best control ###")
            for index_time in eachindex(opt.diagnostics.ρ_best) # sinon pas copy?
                opt.diagnostics.ρ_best[index_time] = copy(opt.model.ρ_t[index_time])
            end
            opt.diagnostics.u_best .= copy(opt.ukp) 
            opt.J_best = opt.J
        end
    end 
    
    opt.iter += 1

    nothing
end

function store_last!(opt::Optim)
    opt.diagnostics.u_last .= copy(opt.ukp)
    for index_time in eachindex(opt.diagnostics.ρ_best) # sinon pas copy?
            opt.diagnostics.ρ_last[index_time] = copy(opt.model.ρ_t[index_time])
    end
    nothing
end 

"""
    diagnostics(opt::Optim; dirname="diagnostics", loss=true, cost=true, staff=true, totvar=false)

TBW
"""
function diagnostics(opt::Optim; dirname="diagnostics", test_case="generic", 
                    loss=true, cost=true, staff=true, totvar=true)
    println("===== diagnostics... =====")
    init_dir(dirname)
    begin
        for (kt, t_snap) in enumerate(opt.saving_times)   #TODO: mettre dans des dossier "/best", "/last"
            println("kt = $kt, t_snap = $t_snap")
            plot_density(opt, kt, t_snap, testcase=test_case, iter2plot="best") 
            plot_density(opt, kt, t_snap, testcase=test_case, iter2plot="last") 
        end

        print("plot mean density...")
        plot_meandensity(opt, iter2plot="best") 
        plot_meandensity(opt, iter2plot="last") 
        println(" done.")

        print("plot control...")
        plot_control(opt, iter2plot="best") 
        plot_control(opt, iter2plot="last") 
        println(" done.")

        print("plot loss history...")
        plot_losshistory(opt, loss=loss, cost=cost, staff=staff, totvar=totvar) 
        println(" done.")

        print("plot residual history...")
        plot_residualhistory(opt) 
        println(" done.")

        print("plot optimality history...")
        plot_optimalityhistory(opt) 
        println(" done.")

        if opt.algo != "FP"
            print("plot descent step history...")
            plot_descentstephistory(opt) 
            println(" done.")
        end
    end
    cd("..") # exit diagnostics directory
    println("===== ...diagnostics done. =====")
end 