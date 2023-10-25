using Plots 
using DelimitedFiles

testcase="1x2" 

loss_GD = readdlm("loss_history_GD")
loss_FP = readdlm("loss_history_FP")
loss_GDFP = readdlm("loss_history_GDFP")

plot(0:length(loss_GD)-1, loss_GD, markershape=:circle, markersize=4, markerstrokewidth=0.5, label="GD", alpha=0.7)
plot!(0:length(loss_FP)-1, loss_FP, markershape=:circle, markersize=4, markerstrokewidth=0.5, label="FP", alpha=0.7)
plot!(0:length(loss_GDFP)-1, loss_GDFP, markershape=:circle, markersize=4, markerstrokewidth=0.5, label="GDFP", alpha=0.7)

# title!("Losses junction "*testcase)
xlabel!("iterations")

savefig("compare_"*testcase*".pdf")

plot!(yscale=:log10, legend=:bottomright)
# title!("Losses junction "*testcase*" (log)")
savefig("compare_"*testcase*"_LOG.pdf")