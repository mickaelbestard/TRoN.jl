export plot_density, plot_control

###############   API   ###############

function plot_losshistory(opt::Optim; loss=true, cost=true, staff=true, totvar=true)

    xlogscale=false
    showmin=true
    linestyle= :solid
    markershape = :circle
    
    label = "loss"
    ylabel = "loss"
    


    DictHistory = Dict()
    if loss DictHistory["loss"]=[copy(opt.diagnostics.losshistory), "steelblue"] end 
    if cost DictHistory["cost"]=[copy(opt.diagnostics.costhistory), "darkorange"] end 
    if staff DictHistory["staff"]=[copy(opt.diagnostics.staffhistory), "purple"] end 
    if totvar DictHistory["TV"]=[copy(opt.diagnostics.totvarhistory), "green"] end 

    ### without coef | not log
    ylogscale=false
    title = "loss history $(opt.algo)"
    plotdicthistory(DictHistory, title=title, xlogscale=xlogscale, ylogscale=ylogscale, showmin=showmin,
    linestyle=linestyle, markershape=markershape)

    ### without coef | log
    ylogscale=true
    title = title * " logscale"
    plotdicthistory(DictHistory, title=title, xlogscale=xlogscale, ylogscale=ylogscale, showmin=showmin,
    linestyle=linestyle, markershape=markershape, legend=:bottomright)


    ### with coef | not log
    θs, θb = opt.regularizations
    if staff DictHistory["staff"][1] .*= θs end 
    if totvar DictHistory["TV"][1] .*= θb end 

    ylogscale=false
    title = "weighted loss history $(opt.algo)"

    plotdicthistory(DictHistory, title=title, xlogscale=xlogscale, ylogscale=ylogscale, showmin=showmin,
    linestyle=linestyle, markershape=markershape)

    ### with coef | log
    ylogscale=true
    title = title * " logscale"

    plotdicthistory(DictHistory, title=title, xlogscale=xlogscale, ylogscale=ylogscale, showmin=showmin,
    linestyle=linestyle, markershape=markershape, legend=:bottomleft)

    nothing
end

function dict2plot_losshistory(D::Dict; loss=true, cost=true, staff=true, totvar=true)

    xlogscale=false
    showmin=true
    linestyle= :solid
    markershape = :circle
    
    label = "loss"
    ylabel = "loss"
    


    DictHistory = Dict()
    if loss DictHistory["loss"]=[copy(D["loss_history"]), "steelblue"] end 
    if cost DictHistory["cost"]=[copy(D["cost_history"]), "darkorange"] end 
    if staff DictHistory["staff"]=[copy(D["staff_history"]), "purple"] end 
    if totvar DictHistory["TV"]=[copy(D["totvar_history"]), "green"] end 

    ### without coef | not log
    ylogscale=false
    title = "loss history $(D["algo"])"
    plotdicthistory(DictHistory, title=title, xlogscale=xlogscale, ylogscale=ylogscale, showmin=showmin,
    linestyle=linestyle, markershape=markershape)

    ### without coef | log
    ylogscale=true
    title = title * " logscale"
    plotdicthistory(DictHistory, title=title, xlogscale=xlogscale, ylogscale=ylogscale, showmin=showmin,
    linestyle=linestyle, markershape=markershape, legend=:bottomright)


    ### with coef | not log
    θs, θb = D["thetaStaff_history"], D["thetaBV_history"]
    if staff DictHistory["staff"][1] .*= θs end 
    if totvar DictHistory["TV"][1] .*= θb end 
    # @show DictHistory["staff"]
    # @show DictHistory["TV"]
    # if staff 
    #     DictHistory["staff"] = DictHistory["staff"] .* θs 
    # end 
    # if totvar 
    #     DictHistory["TV"] = DictHistory["TV"] .* θb 
    # end 

    ylogscale=false
    title = "weighted loss history $(D["algo"])"

    plotdicthistory(DictHistory, title=title, xlogscale=xlogscale, ylogscale=ylogscale, showmin=showmin,
    linestyle=linestyle, markershape=markershape)

    ### with coef | log
    ylogscale=true
    title = title * " logscale"

    plotdicthistory(DictHistory, title=title, xlogscale=xlogscale, ylogscale=ylogscale, showmin=showmin,
    linestyle=linestyle, markershape=markershape, legend=:bottomleft)

    nothing
end

function plot_residualhistory(opt::Optim)

    xlogscale=false
    ylogscale=false
    showmin=true
    linestyle= :solid
    markershape = :circle
    
    label = "residual"
    ylabel = "residual"
    title = "residual history $(opt.algo)"

    history = opt.diagnostics.residualhistory

    start = 2

    plothistory(history, title=title, start=start, label=label, ylabel=ylabel, 
    xlogscale=xlogscale, ylogscale=ylogscale, showmin=showmin,
    linestyle=linestyle, markershape=markershape)

    ### log
    ylogscale=true
    title = title*" logscale"
    plothistory(history, title=title, start=start, label=label, ylabel=ylabel, 
    xlogscale=xlogscale, ylogscale=ylogscale, showmin=showmin,
    linestyle=linestyle, markershape=markershape)

    nothing
end

function dict2plot_residualhistory(D::Dict)

    xlogscale=false
    ylogscale=false
    showmin=true
    linestyle= :solid
    markershape = :circle
    
    label = "residual"
    ylabel = "residual"
    title = "residual history $(D["algo"])"

    history = D["residual_history"]

    start = 2

    plothistory(history, title=title, start=start, label=label, ylabel=ylabel, 
    xlogscale=xlogscale, ylogscale=ylogscale, showmin=showmin,
    linestyle=linestyle, markershape=markershape)

    ### log
    ylogscale=true
    title = title*" logscale"
    plothistory(history, title=title, start=start, label=label, ylabel=ylabel, 
    xlogscale=xlogscale, ylogscale=ylogscale, showmin=showmin,
    linestyle=linestyle, markershape=markershape, legend= :bottomleft)

    nothing
end

function plot_optimalityhistory(opt::Optim)

    xlogscale=false
    ylogscale=false
    showmin=true
    linestyle= :solid
    markershape = :circle
    
    label = L"\Lambda"
    ylabel = "optimality defect"
    title = "optimality history $(opt.algo)"

    history = opt.diagnostics.optimalityhistory

    start = 2

    plothistory(history, title=title, start=start, label=label, ylabel=ylabel, 
    xlogscale=xlogscale, ylogscale=ylogscale, showmin=showmin,
    linestyle=linestyle, markershape=markershape)

    ### log
    ylogscale=true
    title = title*" logscale"
    plothistory(history, title=title, start=start, label=label, ylabel=ylabel, 
    xlogscale=xlogscale, ylogscale=ylogscale, showmin=showmin,
    linestyle=linestyle, markershape=markershape)

    nothing
end

function dict2plot_optimalityhistory(D::Dict)

    xlogscale=false
    ylogscale=false
    showmin=true
    linestyle= :solid
    markershape = :circle
    
    label = L"\Lambda"
    ylabel = "optimality defect"
    title = "optimality history $(D["algo"])"

    history = D["optimality_history"]

    start = 2

    plothistory(history, title=title, start=start, label=label, ylabel=ylabel, 
    xlogscale=xlogscale, ylogscale=ylogscale, showmin=showmin,
    linestyle=linestyle, markershape=markershape)

    ### log
    ylogscale=true
    title = title*" logscale"
    plothistory(history, title=title, start=start, label=label, ylabel=ylabel, 
    xlogscale=xlogscale, ylogscale=ylogscale, showmin=showmin,
    linestyle=linestyle, markershape=markershape, legend= :bottomleft)

    nothing
end

function plot_descentstephistory(opt::Optim)

    xlogscale=false
    ylogscale=false
    showmin=true
    linestyle= :solid
    markershape = :circle
    
    label = L"\delta"
    ylabel = "descent step"
    title = "descent step history $(opt.algo)"

    history = opt.diagnostics.descentstephistory
    start = 2

    plothistory(history, title=title, start=start, label=label, ylabel=ylabel, 
    xlogscale=xlogscale, ylogscale=ylogscale, showmin=showmin,
    linestyle=linestyle, markershape=markershape)

    ### log
    ylogscale=true
    title = title*" logscale"
    plothistory(history, title=title, start=start, label=label, ylabel=ylabel, 
    xlogscale=xlogscale, ylogscale=ylogscale, showmin=showmin,
    linestyle=linestyle, markershape=markershape)

    nothing
end

function dict2plot_descentstephistory(D::Dict)

    xlogscale=false
    ylogscale=false
    showmin=true
    linestyle= :solid
    markershape = :circle
    
    label = L"\delta"
    ylabel = "descent step"
    title = "descent step history $(D["algo"])"

    history = D["descentstep_history"]
    start = 2

    plothistory(history, title=title, start=start, label=label, ylabel=ylabel, 
    xlogscale=xlogscale, ylogscale=ylogscale, showmin=showmin,
    linestyle=linestyle, markershape=markershape)

    ### log
    ylogscale=true
    title = title*" logscale"
    plothistory(history, title=title, start=start, label=label, ylabel=ylabel, 
    xlogscale=xlogscale, ylogscale=ylogscale, showmin=showmin,
    linestyle=linestyle, markershape=markershape)

    nothing
end

"""
    plot_density(opt::Optim, kt, t_snap, path; testcase="generic")

`testcase`` can take the values:
    `generic`   : No assumptions are made about the graph topology. 
                  Each roads are plotted separately. 
    `junction`  : The graph is supposed to be an n×m junction, and n,m are retrieved automatically from the graph. 
                  Roads are plotted according to the graph topology.
    `circle`    : The graph is supposed to be the traffic circle test-case from the paper.
                  Roads are plotted according to the graph topology.
    `threeways` : The graph is supposed to be the threeways test-case from the paper.
                  Roads are plotted according to the graph topology.
"""
function plot_density(opt::Optim, kt, t_snap; testcase="generic", iter2plot="best")
    if testcase == "generic"
        plot_density_generic(opt, kt, t_snap)
    elseif testcase == "junction"
        plot_density_junction(opt, kt, t_snap)
    elseif testcase == "circle"
        plot_density_circle(opt, kt, t_snap)
    elseif testcase == "threeways"
        plot_density_threeways(opt, kt, t_snap, iter2plot=iter2plot)
    else # unknown testcase
        println("Unkown testcase. Switching to generic mode...")
        plot_density_generic(opt::Optim, kt, t_snap)
    end
end

function plot_density(D::Dict, kt, t_snap; testcase="generic", iter2plot="best")
    if testcase == "generic"
        plot_density_generic(D, kt, t_snap, iter2plot=iter2plot)
    elseif testcase == "junction"
        plot_density_junction(D, kt, t_snap, iter2plot=iter2plot)
    elseif testcase == "circle"
        plot_density_circle(D, kt, t_snap, iter2plot=iter2plot)
    elseif testcase == "threeways"
        plot_density_threeways(D, kt, t_snap, iter2plot=iter2plot)
    else # unknown testcase
        println("Unkown testcase. Switching to generic mode...")
        plot_density_generic(D, kt, t_snap, iter2plot=iter2plot)
    end
end

function plot_control(opt::Optim; iter2plot="best")
    ### params 
    title = "Control evolution over time ("*iter2plot*")"
    Nr = opt.model.net.NumberOfRoads
    Nt = opt.model.Nt
    ### control of i-th road is plotted between (xu[i], yu[1]) and (xu[i], yu[2])
    xu = 1:Nr #[1. Nr+2]
    yu = [0 2.]
    steptime=2*yu[2]/opt.model.Tf
    ### TODO: struct ColorParam 
    c1 = colorant"lightblue1"
    c2 = colorant"dodgerblue"
    c3 = colorant"indigo"
    c_step1 = 10
    c_grad = range(c1,c2,length=c_step1+1)
    c_step = 2c_step1
    c_grad = vcat([c_grad[1:end-1] range(c2,c3,length=c_step1)]...) 

    #TODO: mettre dans une autre fonction ?
    #### PLOT
    p = plot(grid= true, ytickfontsize= 7, 
    # left_margin= 2.5mm, right_margin= 3mm, showaxis=:y
    );
    if iter2plot == "best"
        for x in xu # loop over roads
            A = [x yu[1]]; B = [x yu[2]]
            ux = [opt.diagnostics.u_best[t][x] for t in 1:Nt-1]
            plotdashline_u(A, B, Nt-1, ux, c_grad, c_step)
        end
    elseif iter2plot == "last"
        for x in xu # loop over roads
            A = [x yu[1]]; B = [x yu[2]]
            ux = [opt.diagnostics.u_last[t][x] for t in 1:Nt-1]
            plotdashline_u(A, B, Nt-1, ux, c_grad, c_step)
        end
    end
    ### plot colormap
    plotcolormap([1. Nr], [0. 1.], c_grad, c_step, 
                str1=L"u=0", str2=L"u=1", decal=0.15)
    xlabel!("roads") 
    ylabel!(L"time")
    # title!(title)
    xticks!(xu, latexstring.(xu))
    plot!(legend=false)

    f(t; Tf=opt.model.Tf) = (Tf/yu[2]) * t
    # yticks!( 0:40/Nt:2., latexstring.(Int.(round.(f.(0:40/Nt:2.), digits=0))) ) 
    # yticks!( 0:54/Nt:2., latexstring.(round.(f.(0:54/Nt:2.), digits=1))  ) 
    yticks!( 0:steptime:yu[2], latexstring.(round.(f.(0:steptime:yu[2]), digits=1))  ) 
    # annotate!([(0, 2 * 1.1, Plots.text(L"t", 11, :black, :center))])
    # ylims!((-0.1, 2.5))
    # xlims!((-0.05, 6))
    
    savefig(p, "heatmap_control_"*opt.algo*"_"*iter2plot)

    nothing
end

function dict2plot_control(D::Dict; iter2plot="best")
    ### params 
    algo = D["algo"]
    title = "Control evolution over time ("*iter2plot*")"
    Nr = Int(D["NumberOfRoads"][1])
    Nt = Int(D["NumberOfTimeSteps"][1])
    # Tf = Int(D["final_time"][1])
    Tf = D["final_time"][1]
    u_best=D["control_best"]
    u_last=D["control_last"]
    ### control of i-th road is plotted between (xu[i], yu[1]) and (xu[i], yu[2])
    xu = 1:Nr #[1. Nr+2]
    yu = [0 2.]
    steptime=2*yu[2]/Tf
    ### TODO: struct ColorParam 
    c1 = colorant"lightblue1"
    c2 = colorant"dodgerblue"
    c3 = colorant"indigo"
    c_step1 = 10
    c_grad = range(c1,c2,length=c_step1+1)
    c_step = 2c_step1
    c_grad = vcat([c_grad[1:end-1] range(c2,c3,length=c_step1)]...) 

    #TODO: mettre dans une autre fonction ?
    #### PLOT
    p = plot(grid= true, ytickfontsize= 7, 
    # left_margin= 2.5mm, right_margin= 3mm, showaxis=:y
    );
    if iter2plot == "best"
        for x in xu # loop over roads
            A = [x yu[1]]; B = [x yu[2]]
            ux = [u_best[t, x] for t in 1:Nt-1]
            plotdashline_u(A, B, Nt-1, ux, c_grad, c_step)
        end
    elseif iter2plot == "last"
        for x in xu # loop over roads
            A = [x yu[1]]; B = [x yu[2]]
            ux = [u_last[t, x] for t in 1:Nt-1]
            plotdashline_u(A, B, Nt-1, ux, c_grad, c_step)
        end
    end
    ### plot colormap
    plotcolormap([1. Nr], [0. 1.], c_grad, c_step, 
                str1=L"u=0", str2=L"u=1", decal=0.15)
    xlabel!("roads") 
    ylabel!(L"time")
    # title!(title)
    xticks!(xu, latexstring.(xu))
    plot!(legend=false)

    f(t; T=Tf) = (T/yu[2]) * t
    # yticks!( 0:40/Nt:2., latexstring.(Int.(round.(f.(0:40/Nt:2.), digits=0))) ) 
    # yticks!( 0:54/Nt:2., latexstring.(round.(f.(0:54/Nt:2.), digits=1))  ) 
    yticks!( 0:steptime:yu[2], latexstring.(round.(f.(0:steptime:yu[2]), digits=1))  ) 
    # annotate!([(0, 2 * 1.1, Plots.text(L"t", 11, :black, :center))])
    # ylims!((-0.1, 2.5))
    # xlims!((-0.05, 6))
    
    savefig(p, "heatmap_control_"*algo*"_"*iter2plot)

    nothing
end

function plot_meandensity(opt::Optim; iter2plot="best")
    ### params 
    title = "Mean density over time ("*iter2plot*")"
    Nr = opt.model.net.NumberOfRoads
    Nt = opt.model.Nt
    Nc = opt.model.mesh.NumberOfCenters
    ### mean density of i-th road is plotted between (xu[i], yu[1]) and (xu[i], yu[2])
    xr = 1:Nr #[1. Nr+2]
    yr = [0 2.]
    steptime=2*yr[2]/opt.model.Tf
    ### TODO: struct ColorParam 
    c1 = colorant"lime"
    c2 = colorant"blue"
    c3 = colorant"red"
    c_step1 = 10
    c_grad = range(c1,c2,length=c_step1+1)
    c_step = 2c_step1
    c_grad = vcat([c_grad[1:end-1] range(c2,c3,length=c_step1)]...) 

    #TODO: mettre dans une autre fonction ?
    #### PLOT
    p = plot(grid= true, ytickfontsize= 7, 
    # left_margin= 2.5mm, right_margin= 3mm, showaxis=:y
    );
    if iter2plot == "best"
        for r in xr # loop over roads
            A = [r yr[1]]; B = [r yr[2]]
            idx_first = get_global_index_field(opt.model, r, 1)
            idx_last = get_global_index_field(opt.model, r, Nc)
            rhor = [mean(opt.diagnostics.ρ_best[t][idx_first:idx_last]) for t in 1:Nt]
            plotdashline_u(A, B, Nt, rhor, c_grad, c_step)
        end
    elseif iter2plot == "last"
            for r in xr # loop over roads
            A = [r yr[1]]; B = [r yr[2]]
            idx_first = get_global_index_field(opt.model, r, 1)
            idx_last = get_global_index_field(opt.model, r, Nc)
            rhor = [mean(opt.diagnostics.ρ_last[t][idx_first:idx_last]) for t in 1:Nt]
            plotdashline_u(A, B, Nt, rhor, c_grad, c_step)
        end
    end
    
    ### plot colormap
    plotcolormap([1. Nr], [0. 1.], c_grad, c_step, 
                str1=L"\rho_0", str2=L"\rho_{max}", decal=0.15)
    xlabel!("roads") 
    ylabel!(L"time")
    # title!(title)
    xticks!(xr, latexstring.(xr))
    plot!(legend=false)

    f(t; Tf=opt.model.Tf) = (Tf/yr[2]) * t
    # yticks!( 0:40/Nt:2., latexstring.(Int.(round.(f.(0:40/Nt:2.), digits=0))) ) 
    # yticks!( 0:54/Nt:2., latexstring.(round.(f.(0:54/Nt:2.), digits=1))  ) 
    yticks!( 0:steptime:yr[2], latexstring.(round.(f.(0:steptime:yr[2]), digits=1))  ) 
    # annotate!([(0, 2 * 1.1, Plots.text(L"t", 11, :black, :center))])
    # ylims!((-0.1, 2.5))
    # xlims!((-0.05, 6))
    
    savefig(p, "heatmap_meandensity_"*opt.algo*"_"*iter2plot)

    nothing
end

function dict2plot_meandensity(D::Dict; iter2plot="best")
    ### params 
    algo = D["algo"]
    # Tf = Int(D["final_time"][1])
    Tf = D["final_time"][1]
    ρ_best=D["density_best"]
    ρ_last=D["density_last"]
    title = "Mean density over time ("*iter2plot*")"
    Nr = Int(D["NumberOfRoads"][1])
    Nt = Int(D["NumberOfTimeSteps"][1])
    Nc = Int(D["NumberOfCenters"][1])
    ### mean density of i-th road is plotted between (xu[i], yu[1]) and (xu[i], yu[2])
    xr = 1:Nr #[1. Nr+2]
    yr = [0 2.]
    steptime=2*yr[2]/Tf
    ### TODO: struct ColorParam 
    c1 = colorant"lime"
    c2 = colorant"blue"
    c3 = colorant"red"
    c_step1 = 10
    c_grad = range(c1,c2,length=c_step1+1)
    c_step = 2c_step1
    c_grad = vcat([c_grad[1:end-1] range(c2,c3,length=c_step1)]...) 

    #TODO: mettre dans une autre fonction ?
    #### PLOT
    p = plot(grid= true, ytickfontsize= 7, 
    # left_margin= 2.5mm, right_margin= 3mm, showaxis=:y
    );
    if iter2plot == "best"
        for r in xr # loop over roads
            A = [r yr[1]]; B = [r yr[2]]
            idx_first = get_global_index_field(D["adjmat"], r, 1)
            idx_last = get_global_index_field(D["adjmat"], r, :end)
            rhor = [mean(ρ_best[t, idx_first:idx_last]) for t in 1:Nt]
            plotdashline_u(A, B, Nt, rhor, c_grad, c_step)
        end
    elseif iter2plot == "last"
            for r in xr # loop over roads
            A = [r yr[1]]; B = [r yr[2]]
            idx_first = get_global_index_field(D["adjmat"], r, 1)
            idx_last = get_global_index_field(D["adjmat"], r, :end)
            rhor = [mean(ρ_last[t, idx_first:idx_last]) for t in 1:Nt]
            plotdashline_u(A, B, Nt, rhor, c_grad, c_step)
        end
    end
    
    ### plot colormap
    plotcolormap([1. Nr], [0. 1.], c_grad, c_step, 
                str1=L"\rho_0", str2=L"\rho_{max}", decal=0.15)
    xlabel!("roads") 
    ylabel!(L"time")
    # title!(title)
    xticks!(xr, latexstring.(xr))
    plot!(legend=false)

    f(t; Tf=Tf) = (Tf/yr[2]) * t
    # yticks!( 0:40/Nt:2., latexstring.(Int.(round.(f.(0:40/Nt:2.), digits=0))) ) 
    # yticks!( 0:54/Nt:2., latexstring.(round.(f.(0:54/Nt:2.), digits=1))  ) 
    yticks!( 0:steptime:yr[2], latexstring.(round.(f.(0:steptime:yr[2]), digits=1))  ) 
    # annotate!([(0, 2 * 1.1, Plots.text(L"t", 11, :black, :center))])
    # ylims!((-0.1, 2.5))
    # xlims!((-0.05, 6))
    
    savefig(p, "heatmap_meandensity_"*algo*"_"*iter2plot)

    nothing
end

###############   BACKEND   ###############

function plothistory(history; title="loss history GD", 
                        start=1, 
                        ylabel="Jmin", 
                        label="loss",
                        xlogscale=false, 
                        ylogscale=false,
                        linestyle= :solid,
                        markershape = :circle,
                        showmin=true,
                        legend=:topright)

    xx = nothing
    if start == 1 # todo: change name e.g. "start == 0" or "start == :init / :inloop"
        xx = range(start=0, step=1, stop=length(history)-1)
    else 
        xx = range(start=1, step=1, stop=length(history))
    end 

    p=plot(xx, history, label=label, legend=legend, linestyle=linestyle, marker=markershape);

    xlabel!("iterations")
    ylabel!(ylabel) 

    # title!(title)  

    if showmin
        kmin = argmin(history)[1] 
        kmin_plot = xx[kmin] #(start==1) ? kmin : kmin+1
        Hmin = history[kmin_plot]#; println("Hmin = $Hmin")
        str = " at iteration $(kmin_plot)"
        plot!([kmin_plot], [Hmin], seriestype=:scatter, marker=:circle, markersize=4,
            color=:red, label="min : $(round(Hmin, digits=5))"*str);
    end

    xlogscale && plot!(xaxis=:log)
    ylogscale && plot!(yaxis=:log)
    println("at: ", pwd())

    savefig(p, title)
end

function plotdicthistory(dicthistory; 
    title="loss history GD", 
    start=1, 
    ylabel="loss", 
    xlogscale=false, 
    ylogscale=false,
    linestyle= :solid,
    markershape = :circle,
    showmin=true,
    legend=:topright)

    
    xx = nothing
    ratio=1.0

    p=plot(legend=legend, size=ratio .* (600, 420))

    for (key, val) in dicthistory
        
        if start == 1 # todo: change name e.g. "start == 0" or "start == :init / :inloop"
            xx = range(start=0, step=1, stop=length(val[1])-1)
        else 
            xx = range(start=1, step=1, stop=length(val[1]))
        end 

        if minimum(val[1]) == 0 && ylogscale
            val[1] = proj(val[1], a=1e-16, b=maximum(val[1]))
        end
        if key == "loss"
            if ylogscale
                legend = :topright# :bottomleft  
            else 
                legend = :topright 
            end
            plot!(xx, val[1], color=val[2], label=key, linestyle= :solid, legend= legend,
            alpha=0.5, lw=1, marker=:xcross, markersize=4, markeralpha=0.7, markerstrokewidth=0.5)
        else
            plot!(xx, val[1], color=val[2], label=key, linestyle= :dash,# marker= :+, 
            markerstrokewidth=0.1, markersize=4, alpha=0.7)
        end
    end

    xlabel!("iterations")
    ylabel!(ylabel) 

    # title!(title)  
    history=dicthistory["loss"][1]
    if showmin
        kmin = argmin(history)[1] 
        kmin_plot = xx[kmin] #(start==1) ? kmin : kmin-1
        Hmin = history[kmin]#; println("Hmin = $Hmin")
        str = " iter: $(kmin_plot)"
        plot!([kmin_plot], [Hmin], seriestype=:scatter, marker=:circle, markersize=4,
        color=:red, label="min: $(round(Hmin, digits=2))"*str);
    end

    xlogscale && plot!(xaxis=:log)
    ylogscale && plot!(yaxis=:log)
    println("at: ", pwd())

    savefig(p, title * ".pdf")
end

function plot_density_generic(opt::Optim, kt, t_snap; iter2plot="best")
    x = range(0., stop=1.0, length=opt.model.mesh.NumberOfCenters)
    @show index_time = opt.saving_indices[kt]
    for index_road in 1:opt.model.net.NumberOfRoads
        start = get_global_index_field(opt.model, index_road, 1)
        stop  = get_global_index_field(opt.model, index_road, :end)
        if iter2plot=="best"
            rho2plot = opt.diagnostics.ρ_best
        elseif iter2plot=="last"
            rho2plot = opt.diagnostics.ρ_last
        end
        p=plot(x, rho2plot[index_time][start:stop],
        label="density");
        title!("ρ($(round(t_snap, digits=3))) at road $index_road")
        ylims!((-0.01, 1.01))
        namefig = "density"*iter2plot*"_road_$(index_road)_snapshot_$kt.pdf"
        savefig(p, namefig)
    end
    nothing
end

function plot_density_1x1(opt::Optim, kt, t_snap)
    ### parameters 
    filename = "heatmap_1x1_snapshot_$kt.pdf"
    ### vertices are plotted with coordinates (xv[i], yv[i])
    xv = [0 5 10]
    yv = [0 0 0]
    ### edge number rd_num_t[i] is plotted with coordinates (A_t[i], B_t[i])
    A_t = []; B_t = []; rd_num_t = []
    push!(A_t, [0 0], [5 0])
    push!(B_t, [5 0], [10 0])
    push!(rd_num_t, 1, 2)
    index_time = opt.saving_indices[kt]
    density = opt.diagnostics.ρ_best[index_time]
    J = dot(density, opt.redlane)
    xlab = "T = $(round(t_snap, digits=3)),   J(u) = $(round(J, digits=3))"
    plot_heatmap(opt, xv, yv, A_t, B_t, rd_num_t, density,
                filename=filename, xlab=xlab, title="")
    nothing
end

function plot_density_1x2(opt::Optim, kt, t_snap)
    ### parameters 
    filename = "heatmap_1x2_snapshot_$kt"
    #TODOO: struct GraphTopology
    ### vertices are plotted with coordinates (xv[i], yv[i])
    xv = [0 5 10 10] 
    yv = [0 0 5 -5]
    ### edge number rd_num_t[i] is plotted with coordinates (A_t[i], B_t[i])
    A_t = []; B_t = []; rd_num_t = []
    push!(A_t, [0 0], [5 0], [5 0])
    push!(B_t, [5 0], [10 5], [10 -5])
    push!(rd_num_t, 1, 2, 3)
    index_time = opt.saving_indices[kt]
    density = opt.diagnostics.ρ_best[index_time]
    J = dot(density, opt.redlane)
    xlab = "T = $(round(t_snap, digits=3)),   J(u) = $(round(J, digits=3))"
    plot_heatmap(opt, xv, yv, A_t, B_t, rd_num_t, density,
                filename=filename, xlab=xlab, title="")
    nothing
end

function plot_density_2x1(opt::Optim, kt, t_snap)
    ### parameters 
    filename = "heatmap_2x1_snapshot_$kt"
    #TODOO: struct GraphTopology
    ### vertices are plotted with coordinates (xv[i], yv[i])
    xv = [0 0 5 10]
    yv = [5 -5 0 0] 
    ### edge number rd_num_t[i] is plotted with coordinates (A_t[i], B_t[i])
    A_t = []; B_t = []; rd_num_t = []
    push!(A_t, [0 5], [0 -5], [5,0]) 
    push!(B_t, [5 0], [5 0], [10,0])
    push!(rd_num_t, 1, 2, 3)
    index_time = opt.saving_indices[kt]
    density = opt.diagnostics.ρ_best[index_time]
    J = dot(density, opt.redlane)
    xlab = "T = $(round(t_snap, digits=3)),   J(u) = $(round(J, digits=3))"
    plot_heatmap(opt, xv, yv, A_t, B_t, rd_num_t, density,
                filename=filename, xlab=xlab, title="")
    nothing
end

function plot_density_2x2(opt::Optim, kt, t_snap)
    ### parameters 
    filename = "heatmap_2x2_snapshot_$kt"
    #TODOO: struct GraphTopology
    ### vertices are plotted with coordinates (xv[i], yv[i])
    xv = [0 0 5 10 10] 
    yv = [5 -5 0 5 -5]
    ### edge number rd_num_t[i] is plotted with coordinates (A_t[i], B_t[i])
    A_t = []; B_t = []; rd_num_t = []
    push!(A_t, [0 5], [0 -5], [5,0], [5,0]) 
    push!(B_t, [5 0], [5 0], [10,5], [10,-5]) 
    push!(rd_num_t, 1, 2, 3, 4)
    index_time = opt.saving_indices[kt]
    density = opt.diagnostics.ρ_best[index_time]
    J = dot(density, opt.redlane)
    xlab = "T = $(round(t_snap, digits=3)),   J(u) = $(round(J, digits=3))"
    plot_heatmap(opt, xv, yv, A_t, B_t, rd_num_t, density,
                filename=filename, xlab=xlab, title="")
    nothing
end

function plot_density_circle(opt::Optim, kt, t_snap)
    ### parameters 
    filename = "heatmap_circle_snapshot_$(kt)"
    #TODOO: struct GraphTopology
    ### vertices are plotted with coordinates (xv[i], yv[i])
    xv = [0 5 10 10 10 10 15 20]
    yv = [10 10 0 5 15 20 10 10]
    ### edge number rd_num_t[i] is plotted with coordinates (A_t[i], B_t[i])
    A_t = []; B_t = []; rd_num_t = []
    V1 = [9.99 0]; V2 = [10 5]; V3 = [15 10]; V4 = [10 15]; 
    V5 = [10.01 20]; V6 = [20 10]; V7 = [5 10]; V8 = [0 10];
    push!(A_t, V1, V2, V3, V3, V4, V5, V7, V7)
    push!(B_t, V2, V3, V4, V6, V7, V4, V2, V8)
    push!(rd_num_t, 1,2,3,4,5,6,7,8)
    index_time = opt.saving_indices[kt]
    density = opt.diagnostics.ρ_best[index_time]
    J = dot(density, opt.redlane)
    xlab = "T = $(round(t_snap, digits=3)),   J(u) = $(round(J, digits=3))"
    plot_heatmap(opt, xv, yv, A_t, B_t, rd_num_t, density,
                filename=filename, xlab=xlab, title="")
    nothing
end

function plot_density_threeways(opt::Optim, kt, t_snap; iter2plot="best")
    ### parameters 
    filename = "heatmap_threeways_"*iter2plot*"_snapshot_$kt"
    #TODOO: struct GraphTopology
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
    index_time = opt.saving_indices[kt]
    set_control!(opt, opt.diagnostics.u_best)
    ComputePrimal!(opt.model, Constant(0.66))
    if iter2plot=="best"
        density = opt.diagnostics.ρ_best[index_time]
    elseif iter2plot=="last"
        density = opt.diagnostics.ρ_last[index_time]
    end
    J = dot(density, opt.redlane)
    xlab = "T = $(round(t_snap, digits=3)),   J(u) = $(round(J, digits=3))"
    plot_heatmap(opt, xv, yv, A_t, B_t, rd_num_t, density,
                filename=filename, xlab=xlab, title="")
    nothing 
end

function dict2plot_density(D, kt, t_snap; iter2plot="best")
    ### parameters 
    redlane = D["redlane"]

    index_time = Int(D["snapindices"][kt])
    density = D["density_"*iter2plot][index_time, :]

    test_case = D["test_case"]
    filename = "heatmap_"*test_case*"_"*iter2plot*"_snapshot_$kt"

    xv, yv, A_t, B_t, rd_num_t = GraphGeometry[test_case]
    #---
    
    J = dot(density, redlane)
    xlab = "T = $(round(t_snap, digits=3)),   J(u) = $(round(J, digits=3))"

    dict2plot_heatmap(D, xv, yv, A_t, B_t, rd_num_t, density,
                filename=filename, xlab=xlab, title="")
    nothing 
end


function find_junctioncase(opt::Optim)
    nb_ingoing = 0; nb_outgoing = 0
    for v in vertices(opt.model.net.graph)
        if opt.model.net.Junctions[v].Nin * opt.model.net.Junctions[v].Nout != 0
            ### non-leaf junction 
            nb_ingoing = opt.model.net.Junctions[v].Nin
            nb_outgoing = opt.model.net.Junctions[v].Nout
        end
    end
    return nb_ingoing, nb_outgoing
end

function plot_density_junction(opt::Optim, kt, t_snap)
    nb_ingoing, nb_outgoing = find_junctioncase(opt)
    if nb_ingoing == 1
        # 1xm
        if nb_outgoing == 1
            ### 1x1 
            plot_density_1x1(opt, kt, t_snap)
        elseif nb_outgoing == 2
            ### 1x2 
            plot_density_1x2(opt, kt, t_snap)
        end
    elseif nb_ingoing == 2
        # 2xm
        if nb_outgoing == 1
            ### 2x1 
            plot_density_2x1(opt, kt, t_snap)
        elseif nb_outgoing == 2
            ### 2x2 
            plot_density_2x2(opt, kt, t_snap)
        end
    else
        @error "Unknown junction testcase"
    end
end

######################################################################################

# plot a ρ-colored road graph with the corresponding heatmap
function plot_heatmap(opt::Optim, xx, yy, A_t, B_t, rd_num_t, RHO; 
    c1=colorant"lime", c2=colorant"blue", c3=colorant"red", c_step1=15,
    title="T=$(opt.Tf)", filename="", xlab=nothing)
    
    #### COLORS PARAMS
    c_grad = range(c1, c2, length=c_step1+1)
    c_step = 2c_step1
    c_grad = vcat([c_grad[1:end-1] range(c2, c3, length=c_step1)]...)
    
    #### PLOT
    p = plot(axis=([], false));

    # plot edges
    for (A, B, rd_num) in zip(A_t, B_t, rd_num_t)
        plotdashline_ρ(opt, A, B, opt.model.mesh.NumberOfCenters, rd_num, RHO, 
        c_grad, c_step)
    end
    
    # plot nodes
    plot!(xx, yy, seriestype = :scatter, legend=false, markersize=12, color=:gray, alpha=0.5);

    # plot colormap
    plotcolormap(xx, yy, c_grad, c_step)
    !(xlab === nothing) && xlabel!(xlab) 
    title!(title)
    display(p) 
    savefig(p,filename*".pdf") 

    nothing
end

# (from dict) plot a ρ-colored road graph with the corresponding heatmap
function dict2plot_heatmap(D, xx, yy, A_t, B_t, rd_num_t, RHO; 
    c1=colorant"lime", c2=colorant"blue", c3=colorant"red", c_step1=15,
    title="T=$(opt.Tf)", filename="", xlab=nothing)
    
    #### COLORS PARAMS
    c_grad = range(c1, c2, length=c_step1+1)
    c_step = 2c_step1
    c_grad = vcat([c_grad[1:end-1] range(c2, c3, length=c_step1)]...)
    
    #### PLOT
    p = plot(axis=([], false));

    # plot edges
    for (A, B, rd_num) in zip(A_t, B_t, rd_num_t)
        dict2plotdashline_ρ(D, A, B, Int(D["NumberOfCenters"][1]), rd_num, RHO, 
        c_grad, c_step)
    end
    
    # plot nodes
    plot!(xx, yy, seriestype = :scatter, legend=false, markersize=12, color=:gray, alpha=0.5);

    # plot colormap
    plotcolormap(xx, yy, c_grad, c_step)
    !(xlab === nothing) && xlabel!(xlab) 
    title!(title)

    savefig(p,filename*".pdf") 

    nothing
end

# plot a ρ-colored mesh between A and B
function plotdashline_ρ(opt::Optim, A, B, Nc, k, ρ, c_grad, c_step)
    # println("@plotdashline_ρ")
    start = get_global_index_field(opt.model, k, 1)
    stop  = get_global_index_field(opt.model, k, :end)
    ρk = view(ρ, start : stop)
    ρmax = opt.model.net.ρmax_t[k]
    function ρcolor(r) 
        index = max( 1, Int(floor(c_step*r/ρmax)) )
        return c_grad[index] 
    end
    color_tab = ρcolor.(ρk); #println("A = $A, B = $B")
    for j = 1:Nc
        CjL = barycenter(A,B, (j-1)/Nc); #println("CjL = $CjL")
        CjR = barycenter(A,B, j/Nc)    ; #println("CjR = $CjR")
        plotline(CjL, CjR, color_tab[j]) 
    end
end

# (from dict) plot a ρ-colored mesh between A and B
function dict2plotdashline_ρ(D, A, B, Nc, k, ρ, c_grad, c_step)
    # println("@plotdashline_ρ")
    start = get_global_index_field(D["adjmat"], k, 1)
    stop  = get_global_index_field(D["adjmat"], k, :end)
    ρk = view(ρ, start : stop)
    ρmax = D["density_max"][k]
    function ρcolor(r) 
        index = max( 1, Int(floor(c_step*r/ρmax)) )
        return c_grad[index] 
    end
    color_tab = ρcolor.(ρk); #println("A = $A, B = $B")
    for j = 1:Nc
        CjL = barycenter(A,B, (j-1)/Nc); #println("CjL = $CjL")
        CjR = barycenter(A,B, j/Nc)    ; #println("CjR = $CjR")
        plotline(CjL, CjR, color_tab[j]) 
    end
end

# returns the barycenter t, (1-t) between dots A and B 
function barycenter(A,B,t)
    return (1-t)*A .+ t*B
end

# plot the line [A B] with color c
function plotline(A, B, c)
    # println("@plotline()")
    τ = (B[2]-A[2])/(B[1]-A[1]); #println("A = $A, B = $B"); println("τ = $τ")
    x = range(A[1], B[1], length=2)
    y = τ*(x .- A[1]) .+ A[2]
    plot!(x, y, color=c, linewidth=8);

    nothing
end

function plotcolormap(xx, yy, c_grad, c_step; str1=L"0", str2=L"\rho_{max}", 
    decal=0.5, shift=1)
    #shift = 0.7#1.0# 0.1#1.0
    # begin
    #     yy = yy .- 1 ; println("shifting colormap...")#todoo: generic shift
    # end
    xmin  = minimum(xx)
    xmax  = maximum(xx)
    ymin  = minimum(yy) 
    ymax  = maximum(yy) 
    Δy = max(ymax-ymin, 2)/c_step
    ypos(α) = ymin + α*Δy 
    X_ann_left  = 1.1(xmin - decal)
    X_ann_right = 1.1(xmax + decal) +1
    Y_ann_down  = ypos(0) - 0.1*ypos(c_step)#-0.5#0.4#0.1 #decal #(0.1: rho_mean, 0.4: rho_snap)
    Y_ann_up    = ypos(c_step) + 0.1*ypos(c_step) #+0.5#0.4#0.1 #decal

    # xlims!((X_ann_left, X_ann_right+shift))
    # str1== L"0" && ylims!((Y_ann_down, Y_ann_up +0.1 )) #+ shift))
    ylims!((Y_ann_down-0.1, Y_ann_up +0.1 )) #+ shift))

    for j = 1:c_step
        plot!(X_ann_right*ones(2), range(ypos(j-1), ypos(j), length=2), 
              color=c_grad[j], linewidth=10);
    end
    annotate!(X_ann_right,Y_ann_down,text(str1, 10, :black))
    annotate!(X_ann_right,Y_ann_up,text(str2, 10, :black))

    nothing
end


######################################

# plot a u-colored mesh between A and B
function plotdashline_u(A, B, Nt, ux, c_grad, c_step) 
    function ucolor(v) 
        index = max( 1, Int(floor(c_step*v)) )
        return c_grad[index] 
    end
    color_tab = ucolor.(ux)
    for n = 1:Nt#-1
        CnL = barycenter(A,B, (n-1)/Nt) 
        CnR = barycenter(A,B, n/Nt) 
        plot_colorbar(CnL, CnR, color_tab[n]) 
    end

    nothing
end

# plot a color bar as an histogram (to display controls)
function plot_colorbar(A, B, c)
    x = range(A[1], B[1], length=2) #default: length=2
    y = range(A[2], B[2], length=2)
    plot!(x, y, color=c, linewidth=8);

    nothing
end