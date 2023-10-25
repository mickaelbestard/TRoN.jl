using Tron



function test_ComputePrimal(;final_time=0.3, Nc=5)

    network = traffic_circle_v3_graph(1.0)
    model = Model(network, Nc, Tf=final_time)
    init(model)

    initial = Constant(value=0.66)
    ComputePrimal!(model, initial)

    return model
end

function test_ComputeDual(;final_time=0.3, Nc=5, Nmax=2, coef_pen=1e2, coef_reg=1e2, 
    lane=Bool.([1,1,1,0,1,0,0,1]))

    network = traffic_circle_v3_graph(1.0)
    model = Model(network, Nc, Tf=final_time)
    init(model)

    initial = Constant(value=0.66)
    ComputePrimal!(model, initial, computejacobians=true)

    optim = Optim(model, Nmax, coef_pen, coef_reg)
    set_redlane!(optim, lane)

    ComputeDual!(optim)

    return optim
end

function test_Optimizer_circle(;final_time=0.3, Nc=5, Nmax=8, coef_pen=0, coef_reg=0, 
    itermax=3, CV_tol=1e-2, algo="GD", lane=Bool.([1,1,1,0,1,0,0,1]), 
    initstate=Constant(value=0.66), fp_tol=0.0, fp_trigger=10, delay_factor=1.0,
    ismaskcontrol=false, u_init=0.0, lr_init=5e-2, decay=1e-1)


    return test_Optimizer(network=traffic_circle_v3_graph(1.0),
        final_time=final_time, Nc=Nc, Nmax=Nmax, coef_pen=coef_pen, coef_reg=coef_reg,
        itermax=itermax, CV_tol=CV_tol, algo=algo, lane=lane, initstate=initstate,
        fp_tol=fp_tol, fp_trigger=fp_trigger, delay_factor=delay_factor,
        ismaskcontrol=ismaskcontrol, u_init=u_init, lr_init=lr_init, decay=decay)

end

#=function test_Optimizer(;final_time=0.3, Nc=5, Nmax=8, coef_pen=0, coef_reg=0, 
lane=Bool.([1,1,1,0,1,0,0,1]), itermax=3, CV_tol=1e-2, algo="GD", 
initstate=Constant(value=0.66), fp_tol=0.0, fp_trigger=5, delay_factor=1.0,
ismaskcontrol=false, u_init=0.0)

    network = traffic_circle_v3_graph(1.0)
    model = Model(network, Nc, Tf=final_time)
    init(model)

    optim = Optim(model, Nmax, coef_pen, coef_reg)
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

end=#


function compare(;final_time=12, Nc=5, Nmax=2, coef_pen=0, coef_reg=0, 
    itermax=3, CV_tol=1e-2, algo_list=["GD", "FPGD", "GDFP"],
    lane=Bool.([1,1,1,0,1,0,0,1]), 
    initstate=Constant(value=0.66), fp_tol=0.0, fp_trigger=20, delay_factor=1.0,
    ismaskcontrol=false, u_init_list=[0., 0.5, 1.])
    #TODOO: deepcopy of ::Optim --> store in dict (OR ANY OTHER DATA PERSISTANCE)
    # DictOptim = Dict()
    # for algo in algo_list
    #     DictOptim[algo] = []
    # end
    ### assumptions ###
    # u_init : constant
    test_Optimizer(); #trigger compilation
    for u_init in u_init_list
        for algo in algo_list
            println("========== u_init: ", u_init)
            println(">>>>> algo: ", algo)
            optimize = test_Optimizer(final_time=final_time, Nc=Nc, Nmax=Nmax, 
                coef_pen=coef_pen, coef_reg=coef_reg, itermax=itermax, CV_tol=CV_tol, 
                algo=algo, lane=lane, initstate=initstate, fp_tol=fp_tol, 
                fp_trigger=fp_trigger, delay_factor=delay_factor,
                ismaskcontrol=ismaskcontrol, u_init=u_init)
            diagnostics(optimize, test_case="circle", 
                dirname="diagnostics_"*algo*"_u_init_"*string(u_init))
            # push!(DictOptim, optimize)
            optimize = nothing # reset memory ??
            GC.gc(true)
        end
    end 
    # return DictOptim
end

function gridsearch(;final_time=12, Nc=5, Nmax=5, 
    itermax=3, CV_tol=1e-5, algo_list=["GD", "GDFP", "FPGD", "FP"],
    lane=Bool.([1,1,1,0,1,0,0,1]), 
    initstate=Constant(value=0.66), delay_factor=1.0,
    ismaskcontrol=false, u_init_list=[0., 0.5, 1.],
    coef_pen_list=[0, 1e-8, 1e-3, 1e-1],
    coef_reg_list=[0, 1e-8, 1e-3, 1e-1],
    fp_tol_list=[1e-3],
    fp_trigger_list=[10]
    )

    test_Optimizer(); #trigger compilation

    for u_init in u_init_list
        for algo in algo_list
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

                        diagnostics(optimize, test_case="circle", 
                            dirname="diagnostics_"*algo*"_u_init_"*
                            string(u_init)*"_coef_pen_"*string(coef_pen)*
                            "_coef_reg_"*string(coef_reg)*"fp_tol_"*
                            string(fp_tol)*"_fp_trigger_"*string(fp_trigger))
                        # push!(DictOptim, optimize)
                        optimize = nothing # does it reset memory ??
                        GC.gc(true)
                    end 
                end 
            end 
        end
    end 

end

println("load Test_Traffic_Circle") 