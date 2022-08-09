using GasChromatographySimulator
using ForwardDiff
using Plots

# settings
opt = GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-3, "outlet", true)
sys = GasChromatographySimulator.constructor_System(30.0, 0.25e-3, 0.25e-6, "Rxi5MS", "He")
prog =  GasChromatographySimulator.constructor_Program([0.0, 60.0, 600.0, 300.0], [40.0, 40.0, 340.0, 340.0], [401300.0, 401300.0, 401300.0, 401300.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [30.0, 30.0, 30.0, 30.0], [0.0, 0.0, 0.0, 0.0], opt.Tcontrol, sys.L)
sub_C30 = GasChromatographySimulator.Substance("C30",
                                                "638-68-6",
                                                600.477,
                                                37.5137,
                                                181.004,
                                                1e-3,
                                                "Gaida.2021",
                                                5.17619e-5,
                                                0.0,
                                                0.0)
par_C30 = GasChromatographySimulator.Parameters(sys, prog, [sub_C30], opt)

# all calculations for the non-gradient case -> ng=true

# retention time with settings from above:
sol_C30 = GasChromatographySimulator.solve_system_multithreads(par_C30; ng=true)
tR_C30 = sol_C30[1].u[end][1]

function new_thermo_par(Tchar, θchar, ΔCp, par)
    new_sub = GasChromatographySimulator.Substance(par.sub[1].name,
                                                    par.sub[1].CAS,
                                                    Tchar,
                                                    θchar,
                                                    ΔCp,
                                                    par.sub[1].φ₀,
                                                    par.sub[1].ann,
                                                    par.sub[1].Dag,
                                                    par.sub[1].t₀,
                                                    par.sub[1].τ₀)
    new_par = GasChromatographySimulator.Parameters(par.sys, par.prog, [new_sub], par.opt)
    return new_par
end

function tR_thermo4(tp, par)
    # rewrite the solving-function in this form (add parameters L, d, df -> make
    # φ₀=df/d) -> does it work with functions d(x), df(x)?
    # par contains only one solute here?
    function f_tz2(t,p,z)
        Tchar = p[1]
        θchar = p[2]
        ΔCp = p[3]
        φ₀ = p[4]
        f = GasChromatographySimulator.residency(z, t, par.prog.T_itp, par.prog.pin_itp, par.prog.pout_itp, par.sys.L, par.sys.d, par.sys.df, par.sys.gas, ΔCp, Tchar, θchar, φ₀; ng=true)
        return f
    end
    t₀ = par.sub[1].t₀
    zspan = (0.0,par.sys.L)
    prob_tz = ODEProblem(f_tz2, t₀, zspan, tp)
    sol = solve(prob_tz, alg=par.opt.alg, abstol=par.opt.abstol,reltol=par.opt.reltol)
    tR = sol.u[end][1]
    return tR
end
tpar = [590.0, 
        35.0, 
        180.0, 
        1e-3]
tR_thermo4(tpar, par_C30)
tR_thermo4b(tp) = tR_thermo4(tp, par_C30)
tR_thermo4b(tpar)
∂tR4(Tchar, θchar, ΔCp, φ₀) = ForwardDiff.gradient(tR_thermo4b, [Tchar, θchar, ΔCp, φ₀])
∂tR4(590.0, 35.0, 180.0, 1e-3)

# sensitivity analysis
Tchar_range = 550.0:5.0:650.0
θchar_range = 20.0:1.0:40.0
ΔCp_range = 100.0:5.0:200.0
φ₀_range = 0.9e-3:0.01e-3:1.1e-3

tR_var = Array{Float64}(undef, length(Tchar_range), length(θchar_range), length(ΔCp_range), length(φ₀_range))
for i4=1:length(φ₀_range)
    for i3=1:length(ΔCp_range)
        for i2=1:length(θchar_range)
            for i1=1:length(Tchar_range)
                tR_var[i1,i2,i3,i4] = tR_thermo4([Tchar_range[i1], θchar_range[i2], ΔCp_range[i3], φ₀_range[i4]], par_C30)
            end
        end
    end
end
plotly()
plot(Tchar_range, θchar_range, (tR_var[:,:,10,11].-tR_C30).^2, st=:surface)
plot(Tchar_range, ΔCp_range, (tR_var[:,11,:,11].-tR_C30).^2, st=:surface)
plot(Tchar_range, φ₀_range, tR_var[:,11,11,:], st=:surface)

function ressqr(tp)
    return (tR_thermo4(tp, par_C30) - tR_C30)^2
end
∂ressqr(Tchar, θchar, ΔCp, φ₀) = ForwardDiff.gradient(ressqr, [Tchar, θchar, ΔCp, φ₀])
∂ressqr(590.0,35.0, 180.0,1e-3)
∂ressqr(600.477,37.514, 181.0,1e-3)

∂ressqr_var = Array{Float64}(undef, length(Tchar_range), length(θchar_range), length(ΔCp_range), length(φ₀_range))
for i4=1:length(φ₀_range)
    for i3=1:length(ΔCp_range)
        for i2=1:length(θchar_range)
            for i1=1:length(Tchar_range)
                ∂ressqr_var[i1,i2,i3,i4] = ∂ressqr(Tchar_range[i1], θchar_range[i2], ΔCp_range[i3], φ₀_range[i4])[1]
            end
        end
    end
end
plot(Tchar_range, θchar_range, ∂ressqr_var[:,:,10,11], st=:surface)
plot(Tchar_range, ΔCp_range, ∂ressqr_var[:,11,:,11], st=:surface)
plot(Tchar_range, φ₀_range, ∂ressqr_var[:,11,11,:], st=:surface)


function tR_thermo5(tp, par)
    # rewrite the solving-function in this form (add parameters L, d, df -> make
    # φ₀=df/d) -> does it work with functions d(x), df(x)?
    # par contains only one solute here?
    function f_tz2(t,p,z)
        L = p[1]
        a_d = p[2]
        a_df = p[3]
        Tchar = p[4]
        θchar = p[5]
        ΔCp = p[6]
        d(x) = GasChromatographySimulator.gradient(x, [a_d])
        df(x) = GasChromatographySimulator.gradient(x, [a_df])
        φ₀ = a_df/a_d
        f = GasChromatographySimulator.residency(z, t, par.prog.T_itp, par.prog.pin_itp, par.prog.pout_itp, L, d, df, par.sys.gas, ΔCp, Tchar, θchar, φ₀; ng=true)
        return f
    end
    t₀ = par.sub[1].t₀
    zspan = (0.0,par.sys.L)
    prob_tz = ODEProblem(f_tz2, t₀, zspan, tp)
    sol = solve(prob_tz, alg=par.opt.alg, abstol=par.opt.abstol,reltol=par.opt.reltol)
    tR = sol.u[end][1]
    return tR
end

tR_thermo5([30.0, 0.25e-3, 0.25e-6, 600.0, 35.0, 150.0], par_C30)