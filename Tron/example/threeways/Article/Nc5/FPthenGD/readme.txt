opt = test_Optimizer_threeways(final_time=30, Nc=5, algo="FP", fp_tol=1e-3, u_init=1, CV_tol=1e-10, lr_init=0.5, decay=0.01, itermax=500, Nmax=5);

descent_step = opt.scheduler.lr_init / (1 + opt.scheduler.decay*opt.iter)

opt.regularizations .= 1e-8
if opt.res < 1
    opt.regularizations .*= 2
end

if opt.iter > 50 
    optimizer = DictOfOptimizers["GD"]
end


###### junctions infos ######

#1x2 junctions#
alpha = 0.65

#2x1 junctions#
p = 0.7