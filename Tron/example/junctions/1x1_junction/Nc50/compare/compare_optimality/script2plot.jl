using Plots 
using DelimitedFiles

optimality_GD = readdlm("optimality_history_GD")
optimality_FP = readdlm("optimality_history_FP")
optimality_GDFP = readdlm("optimality_history_GDFP")
	
plot(1:length(optimality_GD), optimality_GD, markershape=:circle, markersize=4, markerstrokewidth=0.5, label="GD", alpha=0.7)
plot!(1:length(optimality_FP), optimality_FP, markershape=:circle, markersize=4, markerstrokewidth=0.5, label="FP", alpha=0.7)
plot!(1:length(optimality_GDFP), optimality_GDFP, markershape=:circle, markersize=4, markerstrokewidth=0.5, label="GDFP", alpha=0.7)

# title!("Optimality junction 1x1")
xlabel!("iterations")

savefig("compare_optimality_1x1.pdf")

plot!(yscale=:log10, legend=:bottomright)
# title!("Optimality junction 1x1 (log)")
savefig("compare_optimality_1x1_LOG.pdf")