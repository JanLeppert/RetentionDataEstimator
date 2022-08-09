using GasChromatographySimulator
using ForwardDiff

# -> use the functions of the model in GasChromatographySimulator and redefine the ODE with the expilcit 
# parameters, for which the solutions wil be optimized
# possible parameters:
# - Tchar
# - θchar
# - ΔCp
# - φ₀
# - d
# - L  

#r(x, t, p, par) = GasChromatographySimulator.residency(x, t, par.prog.T_itp, par.prog.Fpin_itp, par.prog.pout_itp, p[1], p[2], p[3], par.col.gas, p[4], p[5], p[6], p[3]/p[2]; ng=par.opt.ng, vis=par.opt.vis, control=par.opt.control)
r(x, t, p, par) = GasChromatographySimulator.residency(x, t, par.prog.T_itp, par.prog.Fpin_itp, par.prog.pout_itp, par.col.L, par.col.d, par.col.df, par.col.gas, p[1], p[2], p[3], p[4]; ng=par.opt.ng, vis=par.opt.vis, control=par.opt.control)

#function rr(x, t, p, par)
#    d(x) = p[2]
#    df(x) = p[3]
#    return GasChromatographySimulator.residency(x, t, par.prog.T_itp, par.prog.Fpin_itp, par.prog.pout_itp, p[1], d(x), df(x), par.col.gas, p[4], p[5], p[6], p[3]/p[2]; ng=par.opt.ng, vis=par.opt.vis, control=par.opt.control)
#end

# example system
theat = [30.0, 60.0, 70.0, 80.0, 90.0, 120.0, 150.0, 180.0]
L = 4.0
d = 0.1e-3
df = 0.1e-6
opt = GasChromatographySimulator.Options(OwrenZen5(), 1e-8, 1e-5, "inlet", true; ng=false, vis="Blumberg", control="Flow")
col = GasChromatographySimulator.Column(L, d, df, "DB5ms", "He")
prog = Array{GasChromatographySimulator.Program}(undef, length(theat))
for i=1:length(theat)
    prog[i] = GasChromatographySimulator.Program([0.0, theat[i]], [40.0, 340.0], [0.8/(60e6), 0.8/(60e6)], [101300.0, 101300.0], [0.0, 0.0], [0.0, 0.0], [col.L, col.L], [0.0, 0.0], opt.Tcontrol, col.L)
end
sub = GasChromatographySimulator.load_solute_database("/Users/janleppert/Documents/GitHub/Misc/plots_Boeker/","oldformat_2022-04-08T16:23:46.009.csv", col.sp, col.gas, ["C10", "C15", "C20", "C25", "C30"], zeros(5), zeros(5))

par = Array{GasChromatographySimulator.Parameters}(undef, length(theat))
for i=1:length(theat)
    par[i] = GasChromatographySimulator.Parameters(col, prog[i], sub, opt)
end

sim = Array{Any}(undef, length(theat))
for i=1:length(theat)
    sim[i] = GasChromatographySimulator.simulate(par[i])
end

retention_times = Array{Float64}(undef, length(sub), length(theat))
for i=1:length(theat)
    retention_times[:,i] = round.(sim[i][1].tR, digits=3) # higher deviation of the retention time by rounding to lesser digits results in quiete different parameter estimations
end

function solving_migration(p, par)
    f_tx(t,p,x) = r(x, t, p, par)
    t₀ = 0.0
    xspan = (0.0, par.col.L)
    prob_tx = ODEProblem(f_tx, t₀, xspan, p)
    solution_tx = solve(prob_tx, alg=par.opt.alg, abstol=par.opt.abstol, reltol=par.opt.reltol)
    return solution_tx
end 

# select the solute
i_s = 1
p0 = [par[1].sub[i_s].Tchar, par[1].sub[i_s].θchar, par[1].sub[i_s].ΔCp, par[1].sub[i_s].φ₀]
test = solving_migration(p0, par[1])

tR_p = Array{Float64}(undef, length(theat))
for i=1:length(theat)
    sol_p = solving_migration(p0, par[i])
    tR_p[i] = sol_p.u[end]
end

using Plots
plotly()
ptR = scatter(theat, retention_times[i_s,:], xlabel="heating time in s", ylabel="tR in s", label="measurement")
scatter!(ptR, theat, tR_p, label="simulation")

function loss_tR(tR, p, par)
    tR_p = Array{Float64}(undef, length(tR))
    for i=1:length(tR)
        sol_p = solving_migration(p, par[i])
        tR_p[i] = sol_p.u[end]
    end
    return sum((tR.-tR_p).^2)
end 

S = loss_tR(retention_times[i_s,:], p0, par)

# what is the idea of this?
p1 = [370.0, 30.0, 100.0, 0.001]
tR_p1 = Array{Float64}(undef, length(theat))
for i=1:length(theat)
    sol_p = solving_migration(p1, par[i])
    tR_p1[i] = sol_p.u[end]
end
scatter!(ptR, theat, tR_p1, label="simulation")

Tchar = 340.0:8.0:420.0
θchar = 20.0:2.0:40.0
ΔCp = 0.0:20.0:200.0

SS = Array{Float64}(undef, length(Tchar), length(θchar), length(ΔCp))
for k=1:length(ΔCp)
    for j=1:length(θchar)
        for i=1:length(Tchar)
            pp = [Tchar[i], θchar[j], ΔCp[k], 0.001]
            SS[i,j,k] = loss_tR(retention_times[i_s,:], pp, par)
        end
    end
end

plot(Tchar, θchar, SS[:,:,5]', st=:surface)
plot(θchar, ΔCp, SS[5,:,:]', st=:surface)
plot(Tchar, ΔCp, SS[:,5,:]', st=:surface)

using Optim

f_opt(x) = loss_tR(retention_times[i_s,:], x, par)

p1 = [sim[2][1].TR[i_s]+273.15, 30.0, 50.0, 0.001]
res_opt_NM = optimize(f_opt, p1, NelderMead())

# add the gradient of f_opt(x) for LBFGS()
res_opt_LBFGS = optimize(f_opt, p1, LBFGS())

p_NM = Optim.minimizer(res_opt_NM)
p_LBFGS = Optim.minimizer(res_opt_LBFGS)

function retention(T, p)
    Tchar = p[1]
    θchar = p[2]
    ΔCp = p[3]
    R = GasChromatographySimulator.R
    lnk = (ΔCp/R + Tchar/θchar)*(Tchar/T-1) + ΔCp/R*log(T/Tchar)
    return lnk
end  

T = 340.0:1.0:420.0
lnk0 = Array{Float64}(undef, length(T)) 
lnk_NM = Array{Float64}(undef, length(T))
lnk_LBFGS = Array{Float64}(undef, length(T))
for i=1:length(T)
    lnk0[i] = retention(T[i], p0)
    lnk_NM[i] = retention(T[i], p_NM)
    lnk_LBFGS[i] = retention(T[i], p_LBFGS)
end 

plnk = plot(T, lnk0, label="original")
plot!(plnk, T, lnk_NM, label="NelderMead")
plot!(plnk, T, lnk_LBFGS, label="LBFGS")

[p0 p_NM p_LBFGS]

S0 = loss_tR(retention_times[i_s,:], p0, par)
S_NM = loss_tR(retention_times[i_s,:], p_NM, par)
S_LBFGS = loss_tR(retention_times[i_s,:], p_LBFGS, par)