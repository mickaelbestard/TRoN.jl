### Explicit time integrators (Runge-Kutta, ...)
### ODE : y'(t) = Flow(y(t), t)

"""
    RK1!(y_new, y_old, Flow, dt)

Generic RK1 step.
"""
function RK1!(y_new, y_old, Flow, dt)
    y_new[:] = y_old .+ dt .* Flow
    nothing
end


########### TODO ###########

# # generic RK2 step
# function RK2!(y_new, y_old, Flow, dt)
#     nothing
# end

# # generic RK3 step
# function RK3!(y_new, y_old, Flow, dt)
#     nothing
# end