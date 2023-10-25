using Tron

function test_ComputePrimal(;final_time=0.3, Nc=5)

    network = graph1x2(1.0)
    model = Model(network, Nc, Tf=final_time)
    init(model)

    initial = Constant(value=0.66)
    ComputePrimal!(model, initial)

    return model
end

function test_ComputeDual(;final_time=0.3, Nc=5, Nmax=2, coef_pen=1e2, coef_reg=1e2, 
    ispenalized=false, isregularized=false, lane=Bool.([1,0,1]))

    network = graph1x2(1.0)
    model = Model(network, Nc, Tf=final_time)
    init(model)

    initial = Constant(value=0.66)
    ComputePrimal!(model, initial, computejacobians=true)

    optim = Optim(model, Nmax, coef_pen, coef_reg, ispenalized, isregularized)
    set_redlane!(optim, lane)
    ComputeDual!(optim)

    return optim
end

function test_Optimizer_1x2(;final_time=0.3, Nc=5, Nmax=2, coef_pen=0, coef_reg=0, 
lane=Bool.([1,0,1]), itermax=3, CV_tol=1e-3, u_init=0.0, delay_factor=1,
algo="GD", initstate=Constant(value=0.66), fp_tol=0.0, fp_trigger=3, ismaskcontrol=false, 
lr_init=5e-2, decay=1e-1)

    return test_Optimizer(network = graph1x2(1.0),
        final_time=final_time, Nc=Nc, Nmax=Nmax, coef_pen=coef_pen, coef_reg=coef_reg,
        itermax=itermax, CV_tol=CV_tol, algo=algo, lane=lane, initstate=initstate,
        fp_tol=fp_tol, fp_trigger=fp_trigger, delay_factor=delay_factor,
        ismaskcontrol=ismaskcontrol, u_init=u_init, lr_init=lr_init, decay=decay)
end

# function test_Optimizer(;final_time=0.3, Nc=5, Nmax=2, coef_pen=0, coef_reg=0, 
# ispenalized=false, isregularized=false, lane=Bool.([1,0,1]), itermax=3, CV_tol=1e-3,
# algo="GD", initstate=Constant(value=0.66), fp_tol=0.0, fp_trigger=3)

#     network = graph1x2(1.0)
#     model = Model(network, Nc, Tf=final_time)
#     init(model)

#     optim = Optim(model, Nmax, coef_pen, coef_reg)
#     set_redlane!(optim, lane)
#     optim.fp_tol = fp_tol
#     optim.trigger = fp_trigger

#     optim(itermax=itermax, CV_tol=CV_tol, algo=algo, initstate=initstate)

#     return optim
# end



println("load Test_1x2")