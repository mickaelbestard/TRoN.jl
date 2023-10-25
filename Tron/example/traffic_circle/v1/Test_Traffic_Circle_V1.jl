using Tron

function test_optim(;Nmax=8, final_time=0.3, Nc=3, itermax=1, optimizer=:grad, momentum=0, plot=false)

    ### graph visualization
    Net = traffic_circle_v1_graph(1., 8) 
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
    @time Test_Traffic_Circle_V1(Tf=final_time, OPTIMIZER=optimizer, MOMENTUM=momentum, itermax=itermax, Nc=Nc, Nmax=Nmax)
end

function test_ODE_primal(;final_time=0.3, Nc=3, plot=false, jac=false)

    ### graph visualization
    Net = traffic_circle_v1_graph(1., 8) 
    if plot 
        g = Net.graph
        edgelabel=["$k" for k in 1:length(edges(g))] 
        gp = gplot(g, edgelabel=edgelabel, nodelabel=vertices(g)) 
        display(gp)
    end
    Opt = OptiNetwork(Net; Nmax=Net.Ne, cfl=0.5, Tf=final_time, Nc=Nc)
    init_Optim(Opt) 
    
    ### run solver
    println("-----------------------")
    println("running ODE_primal()")
    println("final time : $final_time")
    ODE_primal(Opt, jac=jac)

    return Opt
end

function test_ODE_dual(;final_time=0.3, Nc=3, plot=false, 
    lane = Bool.([1,1,1,1,0,0,0,0]))
    
    ### graph visualization
    Net = traffic_circle_v1_graph(1., 8) 
    if plot 
        g = Net.graph
        edgelabel=["$k" for k in 1:length(edges(g))] 
        gp = gplot(g, edgelabel=edgelabel, nodelabel=vertices(g)) 
        display(gp)
    end
    Opt = OptiNetwork(Net; Nmax=8, cfl=0.5, Tf=final_time, Nc=Nc)
    init_Optim(Opt) 
    set_redlane(Opt, lane)
    ### run solver
    println("-----------------------")
    println("running ODE_primal()")
    println("final time : $final_time")
    @time ODE_primal(Opt)
    println("-----------------------")
    println("running ODE_dual()")
    println("final time : $final_time")
    @time ODE_dual(Opt)

    return Opt
end

function diagnostics(;data_path = "./data",
    final_time = 30.,
    optimizer = nothing,
    Net = traffic_circle_v1_graph(1., 8),
    g = Net.graph,
    Nmax = 4,
    Nc = 15,
    algo = "GD")
    
    data_path *= "/lane_11110000/Nmax_$(Nmax)/"

    #TODOO : à mettre en version générique dans Utils.jl
    function extend_lane(lane, Ne=8, Nc=Nc)
        ext_lane = zeros(8*Nc)
        for k in 1:Ne
            for j in 1:Nc 
                ext_lane[(k-1)*Nc + 1 : k*Nc] .= lane[k]
            end
        end
        return ext_lane
    end

    ## vertices
    xr = [0 5 10 10 10 10 15 20]
    yr = [10 10 0 5 15 20 10 10]

    A_t = []; B_t = []; rd_num_t = []
    ## edges #! warning: div/0 in plot for vertical lines
    push!(A_t, [9.99 0], [10 5], [15 10], [10 15], [20 10], [10 15], [0 10], [5 10])
    push!(B_t, [10 5], [15 10], [10 15], [10.01 20], [15 10], [5 10], [5 10], [10 5])
    push!(rd_num_t, 1,2,3,4,5,6,7,8)

    path_t, loss_t, ρs_t, u_t, Λmax_t, norm_Λ_t, first_ctrl_t, last_ctrl_t, total_iter_t, TV_t, Nmax_t = load_data(data_path);

    # gr()
    # timeout = final_time .* [0., 0.2, 0.4, 0.6, 0.8, 1.]  
    timeout = final_time .* collect(0:0.02:1) #[0., 0.2, 0.4, 0.6, 0.8, 1.] 
    #------- CONTROLS -------
    xu = [1. 23] 
    yu = [0 2.];
    #-------------------------
    lane = Bool.([0,0,0,1,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]) 
    ext_lane = extend_lane(lane)
    for k = 1:1#4
        init_dir("dir_$k")
        ### traffic circle v1
        Opt = OptiNetwork(Net; Nmax=Nmax, cfl=0.5, Tf=final_time, Nc=Nc)
        init_Optim(Opt) 
        set_redlane(Opt, lane)
        ######################
        plot_heat_ctrl(Opt, xu, yu, u_t[k], xlab="Roads", filename="ctrl_opt_$(k).png");
        plothistory(loss_t[k], title="loss_$(algo)_$(k)", ylabel="Jmin", label="loss", save=true); 
        plothistory(norm_Λ_t[k], title="||Λ||_$(algo)_$(k)", ylabel="norm(Λ_k)", label="||Λ||_2", save=true);
        plothistory(Λmax_t[k], title="max(|Λ|)_$(algo)_$(k)", ylabel="max(|Λ_k|)", label="max(|Λ|)", save=true);
        plothistory(TV_t[k], title="TotVar_$(k)", ylabel="TV(uk)", save=true);
        plothistory(Nmax_t[k], title="Nmax_$(k)", ylabel="Nmax penalization", save=true);
        time = 0.; nt = 1; n = 1
        for Δt in Opt.Δt_t
            if abs(time - timeout[n]) < Δt
                loss = dot(ext_lane, ρs_t[k][:,nt]) #! loss without penalization / regularization
                filename = "heatmap_traffic_circle_v1_$(algo)_$(n).png"
                xlab = "T=$(round(time,digits=3)),    J(u) = $(round(loss,digits=3))"
                plot_heatmap(Opt, xr, yr, A_t, B_t, rd_num_t, ρs_t[k][:,nt], 
                    filename=filename, xlab=xlab, title="");
                n += 1
            end
            time += Δt
            nt += 1
        end
        if n < length(timeout)+1
            loss = dot(ext_lane, ρs_t[k][:,end])
            filename = "heatmap_traffic_circle_v1_$(algo)_$(n).png"
            xlab = "T=$(timeout[end]),    J(u) = $(round(loss,digits=3))"
            plot_heatmap(Opt, xr, yr, A_t, B_t, rd_num_t, ρs_t[k][:,nt], 
                filename=filename, xlab=xlab, title="");
        end
        init_dir("all_ctrl_$k")
        for it = 1: min(total_iter_t[k], 20)+1
            plot_heat_ctrl(Opt, xu, yu, first_ctrl_t[k][it], xlab="Controlled roads at $(it-1)-th iteration", filename="ctrl_$(k)_iter_$(it-1).png");
            plot_heat_ctrl(Opt, xu, yu, last_ctrl_t[k][end-it+1], xlab="Controlled roads at $(total_iter_t[k]-it+1)-th iteration", filename="ctrl_$(k)_iter_$(total_iter_t[k]-it+1).png");
        end
        cd("..")
        cd("..")
    end
end

println("load Test_Traffic_Circle_V1")