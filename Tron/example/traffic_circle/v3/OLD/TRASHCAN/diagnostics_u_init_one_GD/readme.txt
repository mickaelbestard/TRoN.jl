opt=test_Optimizer(final_time=10, Nc=10, itermax=30, CV_tol=1e-5, u_init=1.0);

if opt.CT<1.0
    descent_step = 10 
else
    norm_Δu = l2_norm(opt.uk .- opt.ukm, opt.model.Δt_t)
    norm_dJ = l2_norm(opt.dJ, opt.model.Δt_t)
    descent_step = max(1000., norm_Δu/(norm_dJ+1e-12))
end

decay_coef = 0.2 
maxiter = 12