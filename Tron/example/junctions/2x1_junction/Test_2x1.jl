using Tron



function test_ComputePrimal(;final_time=0.3, Nc=5)

    network = graph2x1(1.0)
    model = Model(network, Nc, Tf=final_time)
    init(model)

    initial = Constant(value=0.66)
    ComputePrimal!(model, initial)

    return model
end

function test_ComputeDual(;final_time=0.3, Nc=5, Nmax=2, coef_pen=1e2, coef_reg=1e2, 
    ispenalized=false, isregularized=false, lane=Bool.([1,0,1]))

    network = graph2x1(1.0)
    model = Model(network, Nc, Tf=final_time)
    init(model)

    initial = Constant(value=0.66)
    ComputePrimal!(model, initial, computejacobians=true)

    optim = Optim(model, Nmax, coef_pen, coef_reg, ispenalized, isregularized)
    
    set_redlane!(optim, lane)
    ComputeDual!(optim)

    return optim
end

function test_Optimizer_2x1(;final_time=0.3, Nc=5, Nmax=2, coef_pen=0, coef_reg=0, 
lane=Bool.([1,0,1]), itermax=3, CV_tol=1e-3,
algo="GD", initstate=Constant(value=0.66), fp_tol=0.0, fp_trigger=3, delay_factor=1.0,
ismaskcontrol=false, u_init=0.0, lr_init=2e-1, decay=1e-1)

    return test_Optimizer(network = graph2x1(1.0),
        final_time=final_time, Nc=Nc, Nmax=Nmax, coef_pen=coef_pen, coef_reg=coef_reg,
        itermax=itermax, CV_tol=CV_tol, algo=algo, lane=lane, initstate=initstate,
        fp_tol=fp_tol, fp_trigger=fp_trigger, delay_factor=delay_factor,
        ismaskcontrol=ismaskcontrol, u_init=u_init, lr_init=lr_init, decay=decay)

end

# function test_Optimizer(;final_time=0.3, Nc=5, Nmax=2, coef_pen=0, coef_reg=0, 
# ispenalized=false, isregularized=false, lane=Bool.([1,0,1]), itermax=3, CV_tol=1e-3,
# algo="GD", initstate=Constant(value=0.66), fp_tol=0.0, fp_trigger=3)

#     network = graph2x1(1.0)
#     model = Model(network, Nc, Tf=final_time)
#     init(model)

#     optim = Optim(model, Nmax, coef_pen, coef_reg)
#     set_redlane!(optim, lane)
#     optim.fp_tol = fp_tol
#     optim.trigger = fp_trigger

#     optim(itermax=itermax, CV_tol=CV_tol, algo=algo, initstate=initstate)

#     return optim

# end

########### OLD ############


function test_optim(;Nmax=3, final_time=0.3, Nc=3, itermax=1, optimizer=:grad, 
                    momentum=0, plot=false,
                    penalized=false, regularized=false, coef_pen= 1e2, coef_reg= 1e2)
    
    ### graph visualization
    Net = graph2x1(1.)
    if plot 
        g = Net.graph
        edgelabel=["$k" for k in 1:length(edges(g))] 
        gp = gplot(g, edgelabel=edgelabel, nodelabel=vertices(g)) 
        display(gp)
    end

    ### run optimizer
    println("final time : $final_time")
    println("itermax : $itermax")
    println("optimizer : $optimizer")
    println("momentum : $momentum")
    @time Test_2x1(Tf=final_time, OPTIMIZER=optimizer, 
                    MOMENTUM=momentum, itermax=itermax, Nc=Nc,
                    Nmax=Nmax,
                    penalized=penalized, coef_pen=coef_pen,
                    regularized=regularized, coef_reg=coef_reg)
end

function test_ODE_primal(;final_time=0.3, Nc=3, plot=false, jac=false)

    ### graph visualization
    Net = graph2x1(1.) 
    if plot 
        g = Net.graph
        edgelabel=["$k" for k in 1:length(edges(g))] 
        gp = gplot(g, edgelabel=edgelabel, nodelabel=vertices(g)) 
        display(gp)
    end
    Opt = OptiNetwork(Net; Nmax=3, cfl=0.5, Tf=final_time, Nc=Nc)
    init_Optim(Opt) 

    ### run solver
    println("-----------------------")
    println("running ODE_primal()")
    println("final time : $final_time")
    @time ODE_primal(Opt, jac=jac)

    return Opt
end

# function diagnostics(;data_path = "./data",
#     final_time = 0.3,
#     Nmax = 3,
#     Nc = 3,
#     algo = "GD")
    

#     # plotlyjs() #! using PlotlyJS
#     # gr()
#     # pyplot()
    
#     # D = Dict("3"=>"u", "10"=>"û", "17"=>"∇J",
#     #         "20"=>"u\n∇J", "27"=>"um\n∇J")

#     Net = graph2x1(1.)
#     data_path *= "/lane_101/Nmax_$(Nmax)/"

#     #------- DENSITIES -------
#     xr = [0, 0, 5, 10] 
#     yr = [5, -5, 0, 0] 

#     A_t = []; B_t = []; rd_num_t = [] 

#     push!(A_t, [0 5], [0 -5], [5,0]) 
#     push!(B_t, [5 0], [5 0], [10,0])
#     push!(rd_num_t, 1, 2, 3);

#     path_t, loss_t, ρs_t, umin_t, Λmax_t, norm_Λ_t, first_ctrl_t, last_ctrl_t, total_iter_t, TV_t, Nmax_t, res_t, reachmax_t, Λ_t, SΛ_t, U_t, dJ_t, δ_t = load_data(data_path);
#     # for res in res_t
#     #     res[1] = res[2] #? car res[1] arbitraire
#     # end
#     # gr() 
#     timeout = final_time .* collect(0:0.02:1) #[0., 0.2, 0.4, 0.6, 0.8, 1.] 
#     #------- CONTROLS -------
#     xu = [1. 2.]
#     yu = [0 2.];
#     #-------------------------
#     lane = Bool.([1,0,1]) 
#     ext_lane = extend_lane(lane, Net.Ne, Nc)

#     # for k = 1:1#4
#     k=1
#         init_dir("dir_$k")
#         Opt = OptiNetwork(Net; Nmax=Nmax, cfl=0.5, Tf=final_time, Nc=Nc)
#         init_Optim(Opt) 
#         set_redlane(Opt, lane)

#         ### plot optimal control
#         plot_heat_ctrl(Opt, xu, yu, umin_t[k], xlab="Roads", filename="ctrl_opt.png")
#         ### plot loss
#         plothistory(loss_t[k], title="loss_$(algo)", start=0, ylabel="J(u)", label="loss", save=true, ylogscale=false)
#         plothistory(loss_t[k], title="loss_$(algo)_log", start=0, ylabel="J(u)", label="loss", save=true, ylogscale=true)
#         ### plot TV(u)
#         plothistory(TV_t[k], title="TotVar", start=0, ylabel="TV(u)", save=true)
#         ### plot res (= ||Λ||)
#         plothistory(res_t[k], title="res", start=0, ylabel="residual", save=true)
#         ### plot ||Λ|| 
#         plothistory(norm_Λ_t[k], title="||Λ||_$(algo)", ylabel="norm(Λ)", label="||Λ||_2", save=true, ylogscale=false)
#         plothistory(norm_Λ_t[k], title="||Λ||_$(algo)_log", ylabel="norm(Λ)", label="||Λ||_2", save=true, ylogscale=true)
#         ### plot max(Λ)
#         plothistory(Λmax_t[k], title="max(|Λ|)_$(algo)", ylabel="max(|Λ|)", label="max(|Λ|)", save=true, ylogscale=false)
#         plothistory(Λmax_t[k], title="max(|Λ|)_$(algo)_log", ylabel="max(|Λ|)", label="max(|Λ|)", save=true, ylogscale=true)
#         ### plot δ
#         # plothistory(δ_t[k], title="δ", ylabel="grad steps", save=true, ylogscale=false, showmin=false)
#         # plothistory(δ_t[k], title="δ_log", ylabel="grad steps", save=true, ylogscale=true, showmin=false)

#         # plothistory(Nmax_t[k], title="Nmax", ylabel="Nmax penalization", save=true)
#         time = 0.; nt = 1; n = 1
#         for Δt in Opt.Δt_t
#             if abs(time - timeout[n]) < Δt
#                 loss = dot(ext_lane, ρs_t[k][:,nt]) #! loss without penalization / regularization
#                 filename = "heatmap_2x1_$(algo)_$(n).png"
#                 xlab = "T=$(round(time,digits=3)),    J(u) = $(round(loss,digits=3))"
#                 plot_heatmap(Opt, xr, yr, A_t, B_t, rd_num_t, ρs_t[k][:,nt], 
#                     filename=filename, xlab=xlab, title="")
#                 n += 1
#             end
#             time += Δt
#             nt += 1
#         end
#         if n < length(timeout)+1
#             loss = dot(ext_lane, ρs_t[k][:,end])
#             filename = "heatmap_2x1_$(algo)_$(n).png"
#             xlab = "T=$(timeout[end]),    J(u) = $(round(loss,digits=3))"
#             plot_heatmap(Opt, xr, yr, A_t, B_t, rd_num_t, ρs_t[k][:,nt], 
#                 filename=filename, xlab=xlab, title="")
#         end
#         for it in eachindex(Λ_t[k]) 
#             Λ = Λ_t[k][it]
#             SΛ = SΛ_t[k][it]
#             U = U_t[k][it]
#             dJ = dJ_t[k][it]
#             N = size(Λ)[2]
#             times = [0., [sum(Opt.Δt_t[1:nt]) for nt = 1:N-1]...]
#             begin #plot t ↦ Λ3(t); t ↦ [u(t), u(t)-1, ∇J(u)(t)]
#                 Λ3 = Λ[3,:]
#                 U3 = U[3,:]
#                 dJ3 = dJ[3,:]
#                 l = @layout [a b]
#                 p1 = plot(times, Λ3, label="Λ_3", xlabel="time", color=:red, legend= :best)
#                 p3 = plot(times, [U3, U3 .- 1, dJ3], label=["u_3" "u_3 - 1" "dJ(u)_3"], xlabel="time", legend= :best)
#                 p = plot(p1, p3, layout=l, size=(1200,800), ylim=(-1.1,1.1))
#                 savefig(p, "lambda3_$(it)")
#                 # plot!(yaxis=:log)
#                 # savefig(p, "lambda3_$(it)_log")
#             end
#         end

#         # init_dir("all_ctrl_$k")
#         # for it in eachindex(U_t[k])# = 1:min(total_iter_t[k], 20)
#         #     @show it
#         #     plot_heat_ctrl(Opt, xu, yu, first_ctrl_t[k][it], xlab="Controlled roads at $(it-1)-th iteration", filename="ctrl_$(k)_iter_$(it-1).png");
#         #     plot_heat_ctrl(Opt, xu, yu, last_ctrl_t[k][it], xlab="Controlled roads at $(total_iter_t[k]-it+1)-th iteration", filename="ctrl_$(k)_iter_$(total_iter_t[k]-it+1).png");
#         #     # plot_heat_ctrl(Opt, xu, yu, U_t[k][it], xlab="Controlled roads at $(it-1)-th iteration", filename="ctrl_$(k)_iter_$(it-1).png");
#         # end
#         # cd("..")

#         cd("..")
#     # end
# end
println("load Test_2x1")
