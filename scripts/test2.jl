# test of automatic differentiation in ODE solutions with a simple exsample
f(u,p,t) = p*u
p = 1.01
u0 = 1/2
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan,p)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

function sol_of_p(p)
    f(u,p,t) = p*u
    u0 = 1/2
    tspan = (0.0, 1.0)
    prob = ODEProblem(f,u0,tspan,p)
    sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
    return sol.u[end]
end

∂_∂p(p) = ForwardDiff.derivative(sol_of_p, p)

using Plots
plotly()
pp = 0.9:0.01:1.1
plot(pp, sol_of_p.(pp))

plot!(pp, ∂_∂p.(pp))
