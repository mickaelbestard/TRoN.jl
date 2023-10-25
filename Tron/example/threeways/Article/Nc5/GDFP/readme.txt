opt = test_Optimizer_threeways(final_time=30, Nc=5, algo="GDFP", fp_trigger=5, delay_factor=3, fp_tol=1e-3, u_init=1, CV_tol=1e-2, lr_init=1.0, decay=0.01, itermax=500 , Nmax=5);

opt.regularizations .= 1e-8
if opt.res < 1 
    @show opt.regularizations[1] = 1e-4
    @show opt.regularizations[2] = 1e-6
end

descent_step = opt.scheduler.lr_init / (1 + opt.scheduler.decay*opt.iter)

/!\ maskcontrol=true, penser à l'activer manuellement dans Tron.model /!\