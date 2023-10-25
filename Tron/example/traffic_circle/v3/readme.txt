opt = test_Optimizer(final_time=12, 
                    Nc=10, 
                    itermax=200, 
                    algo="GDFPexp", 
                    delay_factor=3, 
                    CV_tol=1e-2, 
                    fp_tol=0);

descent_step = max(1., norm_Δu/(norm_dJ+1e-12))
decay_coef = 0.7 
maxiter = 25

ismaskcontrol = false