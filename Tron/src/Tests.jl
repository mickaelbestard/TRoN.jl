export Test_1x1, Test_1x2, Test_2x1, Test_2x2
export Test_threeways
export Test_Traffic_Circle_V1, Test_Traffic_Circle_V2, Test_Traffic_Circle_V3
export graph1x1, graph1x2, graph2x1, graph2x2
export traffic_circle_graph, threeways_graph


### tests ###



"""
    graph1x1(L_t; rmax_t=[1., 1.], vmax_t=[1., 1.], density_at_leaves_t = [0.66, 0.]  )

    generic writing of the 1x1 test-case
"""
function graph1x1(L_t; rmax_t=[1., 1.], vmax_t=[1., 1.], density_at_leaves_t=[0.66, 0.]  ) 
    
    Nr = 2
    Nv = 3
    Nl = 2

    Net = Network(Nr, Nv, Nl, rmax_t, vmax_t)
    g = Net.graph
    add_edge!(g, 1,2); add_edge!(g, 2,3)

    set_list_leaves(Net, [1,3])
    set_ρ_leaves_t(Net, density_at_leaves_t)
    init_Network(Net)

    # 1x1 junctions
    J1x1 = ones((1,1))
    set_Jmat(Net.Junctions[1], J1x1)
    set_Jmat(Net.Junctions[2], J1x1)
    set_Jmat(Net.Junctions[3], J1x1)

    return Net
end

# function Test_1x1(;L_t=[1., 1.], Nc=3, Tf=1., cfl=0.5, 
#                     LANE=[Bool.([1,1])],
#                     OPTIMIZER=:grad,
#                     MOMENTUM=0.,
#                     itermax=2,
#                     Nmax=2,
#                     penalized=false, coef_pen=1e2,
#                     regularized=false, coef_reg=1e2)

#     init_dir("data")
#     Net = graph1x1(L_t)
#     u_intuitive = zeros(Net.NumberOfRoads)
#     u_intuitive[1] = 1. 
#     println("===== Pipeline : begin =====")
#     tick()
#     Pipeline(Net, LANE, Nmax, cfl, Tf, Nc, 
#         OPTIMIZER=OPTIMIZER, MOMENTUM=MOMENTUM, 
#         itermax=itermax, u_intuitive=u_intuitive,
#         ctrl_lane=true,
#         penalized=penalized, coef_pen=coef_pen,
#         regularized=regularized, coef_reg=coef_reg)
#     tock()
#     println("===== Pipeline : end =====")
# end

# -----------------------------------------------------------------------------

function graph1x2(L_t) # generic writing of the 1x2 test-case
    Nr = 3
    rmax_t = ones(Nr)
    # rmax_t[2] = 0.66
    vmax_t = ones(Nr)
    r_leaves_t = [0.66, 0., 0.] # density at leaves (Dirichlet bc)
    Net = Network(Nr, 4, 3, rmax_t, vmax_t)
    g = Net.graph
    add_edge!(g, 1,2) 
    add_edge!(g, 2,3); add_edge!(g, 2,4)

    set_list_leaves(Net, [1,3,4])
    set_ρ_leaves_t(Net, r_leaves_t)
    init_Network(Net)
    # 1x1 junctions
    J1x1 = ones((1,1))
    set_Jmat(Net.Junctions[1], J1x1)
    set_Jmat(Net.Junctions[3], J1x1)
    set_Jmat(Net.Junctions[4], J1x1)
    # 1x2 junction
    J1x2 = ones((2,1))
    J1x2[1,1]= 0.5 #1. 
    J1x2[2,1]= 0.5 #0. 
    set_Jmat(Net.Junctions[2], J1x2)
    return Net
end

# function Test_1x2(;L=1., Nc=10, Tf=1., cfl=0.5, 
#                     LANE=[Bool.([1,0,1])],
#                     OPTIMIZER=:grad,
#                     MOMENTUM=0.,
#                     itermax=2,
#                     Nmax=3,
#                     penalized=false, coef_pen=1e2,
#                     regularized=false, coef_reg=1e2)
#     init_dir("data")
#     Net = graph1x2(L)
#     u_intuitive = zeros(Net.NumberOfRoads)
#     u_intuitive[1] = 1. 
#     println("===== Pipeline : begin =====")
#     tick()
#     Pipeline(Net, LANE, Nmax, cfl, Tf, Nc, 
#         OPTIMIZER=OPTIMIZER, MOMENTUM=MOMENTUM, 
#         itermax=itermax, u_intuitive=u_intuitive,
#         ctrl_lane=true,
#         penalized=penalized, coef_pen=coef_pen,
#         regularized=regularized, coef_reg=coef_reg)
#     tock()
#     println("===== Pipeline : end =====")
# end

# -----------------------------------------------------------------------------

function graph2x1(L_t) # generic writing of the 2x1 test-case
    Nr = 3
    rmax_t = ones(Nr)
    vmax_t = ones(Nr)
    r_leaves_t = 0.66*[1, 1, 0.] # density at leaves (Dirichlet bc)
    Net = Network(Nr, 4, 3, rmax_t, vmax_t)
    g = Net.graph
    add_edge!(g, 1,3); add_edge!(g, 2,3)
    add_edge!(g, 3,4)

    set_list_leaves(Net, [1,2,4])
    set_ρ_leaves_t(Net, r_leaves_t)
    init_Network(Net)
    # 2x1 junction
    p = 0.7
    J2x1 = ones((1,2))
    set_priorities(Net.Junctions[3], [p]) 
    set_Jmat(Net.Junctions[3], J2x1)
    # 1x1 junctions
    J1x1 = ones((1,1))
    set_Jmat(Net.Junctions[1], J1x1)
    set_Jmat(Net.Junctions[2], J1x1)
    set_Jmat(Net.Junctions[4], J1x1)
    return Net
end

# function Test_2x1(;L=1., Nc=3, Tf=1., cfl=0.5, 
#                     LANE=[Bool.([1,0,1])],
#                     OPTIMIZER=:grad,
#                     MOMENTUM=0.,
#                     itermax=10,
#                     Nmax=3,
#                     penalized=false, coef_pen=1e2,
#                     regularized=false, coef_reg=1e2)

#     init_dir("data")
#     Net = graph2x1(L)
#     u_intuitive = zeros(Net.NumberOfRoads)
#     u_intuitive[1:2] .= 1.
#     println("===== Pipeline : begin =====")
#     tick()
#     Pipeline(Net, LANE, Nmax, cfl, Tf, Nc, 
#         OPTIMIZER=OPTIMIZER, MOMENTUM=MOMENTUM, 
#         itermax=itermax, u_intuitive=u_intuitive,
#         ctrl_lane=true,
#         penalized=penalized, coef_pen=coef_pen,
#         regularized=regularized, coef_reg=coef_reg)
#     tock()
#     println("===== Pipeline : end =====")
# end

# -----------------------------------------------------------------------------

function graph2x2(L_t) # generic writing of the 2x2 test-case
    Nr = 4
    rmax_t = ones(Nr)
    vmax_t = ones(Nr)
    r_leaves_t = [0.66, 0.66, 0., 0.] # density at leaves (Dirichlet bc)
    Net = Network(Nr, 5, 4, rmax_t, vmax_t)
    g = Net.graph
    add_edge!(g, 1,3); add_edge!(g, 3,4)
    add_edge!(g, 2,3); add_edge!(g, 3,5)
    println("Net : $Net")
    set_list_leaves(Net, [1,2,4,5])
    set_ρ_leaves_t(Net, r_leaves_t)
    init_Network(Net)
    # 1x1 junctions (leaves)
    J1x1 = ones((1,1))
    set_Jmat(Net.Junctions[1], J1x1)
    set_Jmat(Net.Junctions[2], J1x1)
    set_Jmat(Net.Junctions[4], J1x1)
    set_Jmat(Net.Junctions[5], J1x1)
    # 2x2 junction
    J2x2 = ones((2,2))
    J2x2[1,1]=0.40; J2x2[1,2]=0.3
    J2x2[2,1]=0.660; J2x2[2,2]=0.7
    set_Jmat(Net.Junctions[3], J2x2)
    return Net
end

# function Test_2x2(;L=1., Nc=10, Tf=2.5, cfl=0.5, 
#                     LANE=[Bool.([1,0,1,0])],
#                     OPTIMIZER=:grad, 
#                     MOMENTUM=0.,
#                     itermax=1,
#                     Nmax=4,
#                     penalized=false, coef_pen=1e2,
#                     regularized=false, coef_reg=1e2)

#     init_dir("data")
#     Net = graph2x2(L)
#     u_intuitive = zeros(Net.NumberOfRoads)
#     u_intuitive[1:2] .= 1.
#     println("===== Pipeline : begin =====")
#     tick()
#     Pipeline(Net, LANE, Nmax, cfl, Tf, Nc, 
#         OPTIMIZER=OPTIMIZER, MOMENTUM=MOMENTUM, 
#         itermax=itermax, u_intuitive=u_intuitive,
#         ctrl_lane=true,
#         penalized=penalized, coef_pen=coef_pen,
#         regularized=regularized, coef_reg=coef_reg)
#     tock()
#     println("===== Pipeline : end =====")
# end

# -----------------------------------------------------------------------------

function traffic_circle_graph(L_t) # generic writing of the traffic circle v2 test-case
    Nr = 8
    rmax_t = ones(Nr)
    vmax_t = ones(Nr)
    r_leaves_t = 0.66 * [1, 1, 0, 0] # density at leaves (Dirichlet bc)
    Net = Network(Nr, 8, 4, rmax_t, vmax_t)
    g = Net.graph
    add_edge!(g, 1,2); add_edge!(g, 2,3); add_edge!(g, 3,4); add_edge!(g, 3,6)
    add_edge!(g, 4,7); add_edge!(g, 5,4); add_edge!(g, 7,2); add_edge!(g, 7,8)
    println("Net : $Net")
    set_list_leaves(Net, [1,5,6,8])
    set_ρ_leaves_t(Net, r_leaves_t)
    init_Network(Net)
    # 1x1 junctions 
    J1x1 = ones(1,1)
    set_Jmat(Net.Junctions[1], J1x1)
    set_Jmat(Net.Junctions[5], J1x1)
    set_Jmat(Net.Junctions[6], J1x1)
    set_Jmat(Net.Junctions[8], J1x1)
    # 1x2 junction
    J1x2 = 0.5*ones(2,1)
    set_Jmat(Net.Junctions[3], J1x2)
    set_Jmat(Net.Junctions[7], J1x2)
    # 2x1 junction
    J2x1 = ones(1,2)
    p = 0.5 ### to give priority on circle: p < 0.5
    set_priorities(Net.Junctions[2], [p, 1-p]) 
    set_Jmat(Net.Junctions[2], J2x1)
    p = 0.5 ### to give priority on circle: p > 0.5
    set_priorities(Net.Junctions[4], [p, 1-p]) 
    set_Jmat(Net.Junctions[4], J2x1)
    return Net
end

# function Test_Traffic_Circle(;L=1., Nc=10, Tf=10., cfl=0.5, 
#                                 LANE=[Bool.([1,1,1,0,1,0,0,1])],
#                                 OPTIMIZER=:grad, 
#                                 MOMENTUM=0.,
#                                 Nmax=8, 
#                                 itermax=1,
#                                 penalized=false, coef_pen=1e2,
#                                 regularized=false, coef_reg=1e2)

#     init_dir("data")
#     Net = traffic_circle_graph(L)
#     u_intuitive = zeros(Net.NumberOfRoads)
#     println("===== Pipeline : begin =====")
#     tick()
#     Pipeline(Net, LANE, Nmax, cfl, Tf, Nc, 
#         OPTIMIZER=OPTIMIZER, MOMENTUM=MOMENTUM, 
#         itermax=itermax, u_intuitive=u_intuitive,
#         ctrl_lane=true,
#         penalized=penalized, coef_pen=coef_pen,
#         regularized=regularized, coef_reg=coef_reg)
#     tock()
#     println("===== Pipeline : end =====")
# end

# -----------------------------------------------------------------------------

function threeways_graph(L_t)
    Nr = 23
    rmax_t = ones(Nr)
    vmax_t = ones(Nr)
    r_leaves_t = [0.66, 0.] # density at leaves (Dirichlet bc)
    Net = Network(Nr, 18, 2, rmax_t, vmax_t)
    g = Net.graph
    ##### graph construction #####
    ## Edges
    # main road
    add_edge!(g, 1,2); add_edge!(g, 2,3); add_edge!(g, 3,4); add_edge!(g, 4,5)
    add_edge!(g, 5,6); add_edge!(g, 6,7); add_edge!(g, 7,8); add_edge!(g, 8,9)
    add_edge!(g, 9,10)
    # upper side
    add_edge!(g, 2,11); add_edge!(g, 11,12); add_edge!(g, 12,13); add_edge!(g, 13,14)
    add_edge!(g, 14,8); add_edge!(g, 12,4); add_edge!(g, 13,6);
    # lower side
    add_edge!(g, 3,15); add_edge!(g, 15,16); add_edge!(g, 16,17); add_edge!(g, 17,18); 
    add_edge!(g, 18,9); add_edge!(g, 16,5); add_edge!(g, 17,7);
    ## Leaves
    set_list_leaves(Net, [1,10])
    set_ρ_leaves_t(Net, r_leaves_t)
    init_Network(Net)
    ## Junctions
    # 1x1 Junctions
    J1x1 = ones((1,1))
    set_Jmat(Net.Junctions[1], J1x1) ; set_Jmat(Net.Junctions[10], J1x1)
    set_Jmat(Net.Junctions[11], J1x1); set_Jmat(Net.Junctions[14], J1x1)
    set_Jmat(Net.Junctions[15], J1x1); set_Jmat(Net.Junctions[18], J1x1)
    # 1x2 junctions
    J1x2 = ones((2,1))
    α = 0.65 
    J1x2[1,1] = α; J1x2[2,1] = 1-α
    set_Jmat(Net.Junctions[2], J1x2) ; set_Jmat(Net.Junctions[3], J1x2)
    set_Jmat(Net.Junctions[12], J1x2); set_Jmat(Net.Junctions[13], J1x2)
    set_Jmat(Net.Junctions[16], J1x2); set_Jmat(Net.Junctions[17], J1x2)
    # 2x1 junction
    p = 0.7
    J2x1 = ones((1,2))
    set_priorities(Net.Junctions[4], [p]); set_priorities(Net.Junctions[5], [p]) 
    set_priorities(Net.Junctions[6], [p]); set_priorities(Net.Junctions[7], [p]) 
    set_priorities(Net.Junctions[8], [p]); set_priorities(Net.Junctions[9], [p]) 
    set_Jmat(Net.Junctions[4], J2x1); set_Jmat(Net.Junctions[5], J2x1) 
    set_Jmat(Net.Junctions[6], J2x1); set_Jmat(Net.Junctions[7], J2x1) 
    set_Jmat(Net.Junctions[8], J2x1); set_Jmat(Net.Junctions[9], J2x1) 
    return Net
end

#LANE=[Bool.([0,1,0,1,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0])]
# function Test_threeways(;L=1., Nc=3, Nr=23, Tf=14., cfl=0.5,
#                         LANE=[Bool.([0,0,0,1,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0])],
#                         OPTIMIZER=:grad, 
#                         MOMENTUM=0.,
#                         itermax=3,
#                         Nmax=10,
#                         penalized=false, coef_pen=1e2,
#                         regularized=false, coef_reg=1e2)

#     init_dir("data")
#     Net = threeways_graph(L)
#     u_intuitive = zeros(Net.NumberOfRoads) 
#     u_intuitive[12:13] .= 1.
#     u_intuitive[18:19] .= 1.
#     println("===== Pipeline : begin =====")
#     tick()
#     Pipeline(Net, LANE, Nmax, cfl, Tf, Nc, 
#         OPTIMIZER=OPTIMIZER, MOMENTUM=MOMENTUM, 
#         itermax=itermax, u_intuitive=u_intuitive,
#         ctrl_lane=false,
#         penalized=penalized, coef_pen=coef_pen,
#         regularized=regularized, coef_reg=coef_reg)
#     tock()
#     println("===== Pipeline : end =====")
# end

# -----------------------------------------------------------------------------
# function traffic_circle_v1_graph(L, Nr) # generic writing of the traffic circle v1 test-case
#     rmax_t = ones(Nr)
#     vmax_t = ones(Nr)
#     r_leaves_t = [0.4, 0.4, 0., 0.4] # density at leaves (Dirichlet bc)
#     Net = Network(Nr, 8, 4, 1, rmax_t, vmax_t)
#     g = Net.graph
#     add_edge!(g, 1,2); add_edge!(g, 2,3); add_edge!(g, 3,4): add_edge!(g, 4,5)
#     add_edge!(g, 4,7); add_edge!(g, 6,3); add_edge!(g, 7,2); add_edge!(g, 8,7)
#     println("Net : $Net")
#     set_list_leaves(Net, [1,5,6,8])
#     set_ρ_leaves_t(Net, r_leaves_t)
#     init_Network(Net, L)
#     # 1x1 junctions 
#     J1x1 = ones(1,1)
#     set_Jmat(Net.Junctions[1], J1x1)
#     set_Jmat(Net.Junctions[5], J1x1)
#     set_Jmat(Net.Junctions[6], J1x1)
#     set_Jmat(Net.Junctions[8], J1x1)
#     # 1x2 junction
#     J1x2 = 0.5*ones(2,1)
#     set_Jmat(Net.Junctions[4], J1x2)
#     # 2x1 junction
#     p = 0.5
#     J2x1 = ones(1,2)
#     set_priorities(Net.Junctions[2], [p]) 
#     set_Jmat(Net.Junctions[2], J2x1)
#     set_priorities(Net.Junctions[3], [p]) 
#     set_Jmat(Net.Junctions[3], J2x1)
#     set_priorities(Net.Junctions[7], [p]) 
#     set_Jmat(Net.Junctions[7], J2x1)
#     return Net
# end

# function Test_Traffic_Circle_V1(;Nr=8, L=1., Nc=10, Tf=10., cfl=0.5, 
#         LANE=[Bool.([1,1,1,1,0,0,0,0])],
#         OPTIMIZER=:grad, 
#         MOMENTUM=0.,
#         Nmax=8, 
#         itermax=1)

#     init_dir("data")
#     Net = traffic_circle_v1_graph(L, Nr)
#     u_intuitive = zeros(Net.NumberOfRoads)
#     # u_intuitive[1:2] .= 1.
#     println("===== Pipeline : begin =====")
#     tick()
#     Pipeline(Net, LANE, Nmax, cfl, Tf, Nc, 
#     OPTIMIZER=OPTIMIZER, MOMENTUM=MOMENTUM, 
#     itermax=itermax, u_intuitive=u_intuitive,
#     ctrl_lane=false)
#     tock()
#     println("===== Pipeline : end =====")
# end

# # -----------------------------------------------------------------------------
# function traffic_circle_v2_graph(L, Nr) # generic writing of the traffic circle v2 test-case
#     rmax_t = ones(Nr)
#     vmax_t = ones(Nr)
#     r_leaves_t = [0.4, 0., 0.4, 0.] # density at leaves (Dirichlet bc)
#     Net = Network(Nr, 8, 4, 1, rmax_t, vmax_t)
#     g = Net.graph
#     add_edge!(g, 1,2); add_edge!(g, 2,3); add_edge!(g, 3,4): add_edge!(g, 4,5)
#     add_edge!(g, 4,7); add_edge!(g, 6,3); add_edge!(g, 7,2); add_edge!(g, 7,8)
#     println("Net : $Net")
#     set_list_leaves(Net, [1,5,6,8])
#     set_ρ_leaves_t(Net, r_leaves_t)
#     init_Network(Net, L)
#     # 1x1 junctions 
#     J1x1 = ones(1,1)
#     set_Jmat(Net.Junctions[1], J1x1)
#     set_Jmat(Net.Junctions[5], J1x1)
#     set_Jmat(Net.Junctions[6], J1x1)
#     set_Jmat(Net.Junctions[8], J1x1)
#     # 1x2 junction
#     J1x2 = 0.5*ones(2,1)
#     set_Jmat(Net.Junctions[4], J1x2)
#     set_Jmat(Net.Junctions[7], J1x2)
#     # 2x1 junction
#     p = 0.5
#     J2x1 = ones(1,2)
#     set_priorities(Net.Junctions[2], [p]) 
#     set_Jmat(Net.Junctions[2], J2x1)
#     set_priorities(Net.Junctions[3], [p]) 
#     set_Jmat(Net.Junctions[3], J2x1)
#     return Net
# end

# function Test_Traffic_Circle_V2(;Nr=8, L=1., Nc=10, Tf=10., cfl=0.5, 
#         LANE=[Bool.([1,1,1,1,0,0,0,0])],
#         OPTIMIZER=:grad, 
#         MOMENTUM=0.,
#         Nmax=8, 
#         itermax=1,
#         penalized=penalized, coef_pen=coef_pen,
#         regularized=regularized, coef_reg=coef_reg)

#     init_dir("data")
#     Net = traffic_circle_v2_graph(L, Nr)
#     u_intuitive = zeros(Net.NumberOfRoads)
#     # u_intuitive[1:2] .= 1.
#     println("===== Pipeline : begin =====")
#     tick()
#     Pipeline(Net, LANE, Nmax, cfl, Tf, Nc, 
#         OPTIMIZER=OPTIMIZER, MOMENTUM=MOMENTUM, 
#         itermax=itermax, u_intuitive=u_intuitive,
#         ctrl_lane=false,
#         penalized=penalized, coef_pen=coef_pen,
#         regularized=regularized, coef_reg=coef_reg)
#     tock()
#     println("===== Pipeline : end =====")
# end

# -----------------------------------------------------------------------------


println("load Tests.jl")