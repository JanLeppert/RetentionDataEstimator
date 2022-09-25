using ForwardDiff
using GasChromatographySimulator
using Plots
using Optimization
using OptimizationOptimJL
using OptimizationBBO
using Dierckx
using BenchmarkTools

include("/Users/janleppert/Documents/GitHub/ThermodynamicDataEstimator/src/ThermodynamicDataEstimator.jl")

# test system (put this in a function)
# Definitions
db_path = "/Users/janleppert/Documents/GitHub/GasChromatographySimulator/data"
db_file = "Database_test.csv"
solutes = ["Octane", "2-Octanone", "2-Octanol", "1-Octanol"]
Tinit = 40.0 # start temperature in °C
tinit = 1.0 # hold time of start temperature in min
Tend = 280.0 # end temperature in °C
tend = 2.0 # hold time of the end temperature in min
heating_rates = [1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0] # °C/min 
pinit = 150000.0
pend = 250000.0
TPs = Array{Array{Float64,1}}(undef, length(heating_rates))
PPs = Array{Array{Float64,1}}(undef, length(heating_rates))
for i=1:length(heating_rates)
	TPs[i] = [Tinit, tinit, heating_rates[i], Tend, tend]
	pressure_rate = (pend-pinit)/(Tend-Tinit)*heating_rates[i]
	PPs[i] = [pinit, tinit, pressure_rate, pend, tend]
end
L = 30.0
d = 0.25e-3
df = 0.25e-6
sp = "SLB5ms"
gas = "He"

# simulate test data
tR_meas, tR_meas_rand, par_meas = ThermodynamicDataEstimator.sim_test_chrom(L, d, df, sp, gas, TPs, PPs, solutes, db_path, db_file)


## here the estimation of the parameters begins
# 0. Data of the used GC system
#       - for real systems, use measured parameters, if possible
L_def = L
d_def = d
df_def = df
heating_rates_def = heating_rates
TPs_def = TPs
PPs_def = PPs
par_def = par_meas

prog_def = Array{GasChromatographySimulator.Program}(undef, length(TPs_def))
for i=1:length(TPs_def)
    prog_def[i] = GasChromatographySimulator.Program(TPs_def[i], PPs_def[i], L_def; pout="vacuum", time_unit="min")
end

# 1. Start values of the parameters
Tchar_est, θchar_est, ΔCp_est = ThermodynamicDataEstimator.estimate_start_parameter(tR_meas, TPs_def, PPs_def, L_def, d_def, gas; pout="vacuum", time_unit="min", control="Pressure")

# true values of the Kcentric parameters
Tchar = Array{Float64}(undef, length(solutes))
θchar = Array{Float64}(undef, length(solutes))
ΔCp = Array{Float64}(undef, length(solutes))
for i=1:length(solutes)
    Tchar[i] = par_meas[1].sub[i].Tchar
    θchar[i] = par_meas[1].sub[i].θchar
    ΔCp[i] = par_meas[1].sub[i].ΔCp
end

# 2. Optimization
# optimize every solute separatly
function optimize_Kcentric_single(tR, L, d, gas, prog, opt, Tchar_e, θchar_e, ΔCp_e, method, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp)
	optf = OptimizationFunction(ThermodynamicDataEstimator.opt_Kcentric, Optimization.AutoForwardDiff())
	
    opt_sol = Array{Any}(undef, size(tR)[2])
	for i=1:size(tR)[2]
		p = [tR[:,i], L, d, prog, opt, gas]
		x0 = [Tchar_e[i], θchar_e[i], ΔCp_e[i]]
		lb = [lb_Tchar[i], lb_θchar[i], lb_ΔCp[i]]
		ub = [ub_Tchar[i], ub_θchar[i], ub_ΔCp[i]]
		if method == NelderMead() || method == NewtonTrustRegion() || Symbol(method) == Symbol(Newton())
			prob = OptimizationProblem(optf, x0, p)
		else
			prob = OptimizationProblem(optf, x0, p, lb=lb, ub=ub)
		end
		opt_sol[i] = solve(prob, method)
	end
	return opt_sol
end

opt = GasChromatographySimulator.Options(ng=true)

methods = [NelderMead(), NewtonTrustRegion(), Newton()] # add more methods
# exact retention times
sol_single = Array{Any}(undef, length(methods))
for i=1:length(methods)
    sol_single[i] = optimize_Kcentric_single(tR_meas, L_def, d_def, "He", prog_def, opt, Tchar_est, θchar_est, ΔCp_est, methods[i], Tchar_est*0.9, θchar_est*0.8, ΔCp_est*0.5, Tchar_est*1.1, θchar_est*1.2, ΔCp_est*1.5)
end
# shifted retention times
sol_single_rand = Array{Any}(undef, length(methods))
for i=1:length(methods)
    sol_single_rand[i] = optimize_Kcentric_single(tR_meas_rand, L_def, d_def, "He", prog_def, opt, Tchar_est, θchar_est, ΔCp_est, methods[i], Tchar_est*0.9, θchar_est*0.8, ΔCp_est*0.5, Tchar_est*1.1, θchar_est*1.2, ΔCp_est*1.5)
end

# benchmark time measurement 
b_single = Array{Any}(undef, length(methods))
b_single_rand = Array{Any}(undef, length(methods))
#for i=1:length(methods)
    b_single[1] = @benchmark optimize_Kcentric_single(tR_meas, L_def, d_def, "He", prog_def, opt, Tchar_est, θchar_est, ΔCp_est, methods[1], Tchar_est*0.9, θchar_est*0.8, ΔCp_est*0.5, Tchar_est*1.1, θchar_est*1.2, ΔCp_est*1.5)
    b_single[2] = @benchmark optimize_Kcentric_single(tR_meas, L_def, d_def, "He", prog_def, opt, Tchar_est, θchar_est, ΔCp_est, methods[2], Tchar_est*0.9, θchar_est*0.8, ΔCp_est*0.5, Tchar_est*1.1, θchar_est*1.2, ΔCp_est*1.5)
    b_single[3] = @benchmark optimize_Kcentric_single(tR_meas, L_def, d_def, "He", prog_def, opt, Tchar_est, θchar_est, ΔCp_est, methods[3], Tchar_est*0.9, θchar_est*0.8, ΔCp_est*0.5, Tchar_est*1.1, θchar_est*1.2, ΔCp_est*1.5)
    b_single_rand[1] = @benchmark optimize_Kcentric_single(tR_meas_rand, L_def, d_def, "He", prog_def, opt, Tchar_est, θchar_est, ΔCp_est, methods[1], Tchar_est*0.9, θchar_est*0.8, ΔCp_est*0.5, Tchar_est*1.1, θchar_est*1.2, ΔCp_est*1.5)
    b_single_rand[2] = @benchmark optimize_Kcentric_single(tR_meas_rand, L_def, d_def, "He", prog_def, opt, Tchar_est, θchar_est, ΔCp_est, methods[2], Tchar_est*0.9, θchar_est*0.8, ΔCp_est*0.5, Tchar_est*1.1, θchar_est*1.2, ΔCp_est*1.5)
    b_single_rand[3] = @benchmark optimize_Kcentric_single(tR_meas_rand, L_def, d_def, "He", prog_def, opt, Tchar_est, θchar_est, ΔCp_est, methods[3], Tchar_est*0.9, θchar_est*0.8, ΔCp_est*0.5, Tchar_est*1.1, θchar_est*1.2, ΔCp_est*1.5)
#end

# save/export the results and benchmark times to .csv 