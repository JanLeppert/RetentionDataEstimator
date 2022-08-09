### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ bdaf89aa-42cd-4ed4-b009-a7820e233cc9
begin
	import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
	include("/Users/janleppert/Documents/GitHub/ThermodynamicDataEstimator/src/ThermodynamicDataEstimator.jl")
	using ForwardDiff
	using GasChromatographySimulator
	using Optim
	using Plots
	using Roots
	using IntervalRootFinding
	using StaticArrays
	using LeastSquaresOptim
	using BlackBoxOptim
	using PlutoUI
	TableOfContents()
end

# ╔═╡ 2c5e497e-09f2-4ffe-8df5-c455fd275771
md"""
# Test of Estimation of K-centric Parameters
"""

# ╔═╡ f3ab4008-8851-4ba2-87bc-eb54ebbc0d7d
md"""
## Construct a test system
"""

# ╔═╡ 02ec232f-83f8-4b60-afed-904cdc08e106
# Temperature Program
TP = [40.0, 1.0, 5.0, 200.0, 0.0, 10.0, 280.0, 2.0, 20.0, 320.0, 2.0]

# ╔═╡ fdecf17b-0b87-4df0-9bb5-a946505ce290
# Pressure Program
PP = [200000.0, 10.0, 3000.0, 300000.0, 5.0]

# ╔═╡ d9473aa5-40cd-47aa-ad5a-76c9ebfc1ade
db_path = "/Users/janleppert/Documents/GitHub/GasChromatographySimulator/data"

# ╔═╡ 4289f252-d7ba-4a21-aea3-1555e3d22d22
db_file = "Database_test.csv"

# ╔═╡ dab67cfa-dd79-44d4-bd78-327dafb550ed
solutes = ["Octane", "2-Octanone", "1-Octanol", "2-Octanol"]

# ╔═╡ 8bc58767-2d14-4317-b94b-8c47e2564f30
# function to create a conventional GC-system
function conventional_GC(L, d, df, sp, gas, TP, PP, solutes, db_path, db_file)
	opt = GasChromatographySimulator.Options()
	col = GasChromatographySimulator.Column(L, d, df, sp, gas)
	prog = GasChromatographySimulator.Program(TP, PP, L; pout="vacuum", time_unit="min")
	sub = GasChromatographySimulator.load_solute_database(db_path, db_file, sp, gas, solutes, zeros(length(solutes)), zeros(length(solutes)))
	par = GasChromatographySimulator.Parameters(col, prog, sub, opt)
	return par
end

# ╔═╡ c2a06f83-0c27-4841-ba3e-fd8eb34565a2
par_test = conventional_GC(30.0, 0.25e-3, 0.25e-6, "SLB5ms", "He", TP, PP, solutes, db_path, db_file)

# ╔═╡ 33d48cb9-a06c-4669-9df3-a93de198e8a9
GasChromatographySimulator.plot_temperature(par_test)

# ╔═╡ 6039f94d-baf2-4757-9779-29791b1addc9
plot(GasChromatographySimulator.plot_pressure(par_test), GasChromatographySimulator.plot_flow(par_test))

# ╔═╡ 480a7a33-d13c-40ce-83ee-28ff647e7a70
md"""
## Simulate the test system
"""

# ╔═╡ e8368d40-faae-4a4c-b88a-985da383a4c0
md"""
### Using ODEsys
"""

# ╔═╡ 5f9870a7-28c8-49a8-9ba6-77b94df6c908
pl_test, sol_test = GasChromatographySimulator.simulate(par_test)

# ╔═╡ aac1d065-18a5-4109-8dc2-510b3a6cab0c
GasChromatographySimulator.plot_chromatogram(pl_test, [0.0, 2000.0])[1]

# ╔═╡ d95f2b44-f7b7-4d6c-a9c8-7b2d5b462b37
md"""
### Only migration
"""

# ╔═╡ c22086bf-17c5-4608-8e44-26a30f0689d2
# simulate, using only the migration function
begin
	sol_tz = Array{Any}(undef, length(solutes))
	tR = Array{Float64}(undef, length(solutes)) 
	for i=1:length(solutes)
		sol_tz[i] = GasChromatographySimulator.solving_migration(par_test.col, par_test.prog, par_test.sub[i], par_test.opt)
		tR[i] = sol_tz[i].u[end]
	end
end

# ╔═╡ 1cd045bd-4d72-4b4d-94be-ceb452f7e989
tR

# ╔═╡ c5b85477-ed0a-42b7-bd4e-aa1282871b19
ΔtR = pl_test.tR .- tR

# ╔═╡ 4fcc7a44-7899-45ac-a4d6-b1e53bec411d
md"""
### Note

The difference between the two simulation modes (ODEsystem and migration only) result in nearly the same retention times. Differences are neglectable (max difference = $(round(maximum(abs.(ΔtR)); sigdigits=3))s).
"""

# ╔═╡ 708d03c2-546e-4f89-be85-f9d57e656a01
md"""
## Create measurements
"""

# ╔═╡ 94424a15-c0fb-4f0c-be1c-6146585f1e25
Tinit = 40.0 # start temperature in °C

# ╔═╡ 4919d447-709e-4261-acb2-e713937f86fc
tinit = 1.0 # hold time of start temperature in min

# ╔═╡ 00603cca-e1bd-4db7-af0f-a7661edbc109
Tend = 280.0 # end temperature in °C

# ╔═╡ 34dd3b26-4ef8-47da-837a-a6daee7356d8
tend = 2.0 # hold time of the end temperature in min

# ╔═╡ 9b9f492b-d4aa-4ac1-9906-51b3f29ce17d
heating_rates = [1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0] # °C/min 

# ╔═╡ 9ccf93d6-0426-4d33-86f2-21537811bfcb
pinit = 150000.0

# ╔═╡ 98cb52b0-65a4-4fde-836b-e8837363c00e
pend = 250000.0

# ╔═╡ ed18e731-67ba-433e-a802-a258f097bca4
begin
	TPs = Array{Array{Float64,1}}(undef, length(heating_rates))
	PPs = Array{Array{Float64,1}}(undef, length(heating_rates))
	for i=1:length(heating_rates)
		TPs[i] = [Tinit, tinit, heating_rates[i], Tend, tend]
		pressure_rate = (pend-pinit)/(Tend-Tinit)*heating_rates[i]
		PPs[i] = [pinit, tinit, pressure_rate, pend, tend]
	end
end

# ╔═╡ 4d144797-380b-45e0-b412-63b08bd21762
TPs

# ╔═╡ f48310e7-8ffc-43cd-be23-169c418ba80d
PPs

# ╔═╡ 906da7f5-ba06-4ecd-9850-084a938cafc4
begin
	par_meas = Array{GasChromatographySimulator.Parameters}(undef, length(heating_rates))
	prog_def = Array{GasChromatographySimulator.Program}(undef, length(heating_rates))
	for i=1:length(heating_rates)
		par_meas[i] = conventional_GC(30.0, 0.25e-3, 0.25e-6, "SLB5ms", "He", TPs[i], PPs[i], solutes, db_path, db_file)
		prog_def[i] = par_meas[i].prog
	end
	col_def = par_meas[1].col
	sub_def = par_meas[1].sub
	opt_def = par_meas[1].opt
end

# ╔═╡ 915a8459-2265-4a7d-9113-de562dfc02a5
ii = 7

# ╔═╡ 322d274a-bfc7-46af-999a-6d3076f905c3
plot(GasChromatographySimulator.plot_temperature(par_meas[ii]), GasChromatographySimulator.plot_flow(par_meas[ii]))

# ╔═╡ 2446c0d3-8105-4644-a67d-ec39ec4dc1fc
begin
	pl_meas = Array{DataFrame}(undef, length(heating_rates))
	tR_meas = Array{Float64}(undef, length(heating_rates), length(solutes))
	tR_meas_randn = Array{Float64}(undef, length(heating_rates), length(solutes))
	measurements = DataFrame(heating_rate=heating_rates)
	measurements_randn = DataFrame(heating_rate=heating_rates)
	for i=1:length(heating_rates)
		pl_meas[i] = GasChromatographySimulator.simulate(par_meas[i])[1]
		for j=1:length(solutes)
			jj = findfirst(pl_meas[i].Name.==solutes[j])
			tR_meas[i,j] = pl_meas[i].tR[jj]
			tR_meas_randn[i,j] = pl_meas[i].tR[jj]*(1+randn()*0.005) # add here a random ± time 
		end
	end
	for j=1:length(solutes)
		measurements[!, Symbol(string("tR_", solutes[j]))] = tR_meas[:,j]
		measurements_randn[!, Symbol(string("tR_", solutes[j]))] = tR_meas_randn[:,j]
	end
	measurements, measurements_randn
end

# ╔═╡ 5b136b18-5a0f-4cff-8a1c-585b390bab0c
188.363 + randn()*188.363*0.005

# ╔═╡ f1c1d843-ef16-4212-81c1-56ede0301c54
188.363*(1+randn()*0.01)

# ╔═╡ 5ce58c17-ec49-4fa1-abed-961cb1380d6a
pl_meas

# ╔═╡ 599cf329-ed37-453c-a190-b909847c1ec9
md"""
## Define functions for optimization
"""

# ╔═╡ 7011dff6-2e07-4d03-b21f-94fd4359adc6
r_Kcentric(x, t, p, par_def) = GasChromatographySimulator.residency(x, t, par_def.prog.T_itp, par_def.prog.Fpin_itp, par_def.prog.pout_itp, par_def.col.L, par_def.col.d, par_def.col.df, par_def.col.gas, p[1], p[2], p[3], par_def.col.df/par_def.col.d; ng=par_def.opt.ng, vis=par_def.opt.vis, control=par_def.opt.control) # p[1] -> Tchar, p[2] -> θchar, p[3] -> ΔCp

# ╔═╡ 9fbd5c15-1bb5-4344-b6a4-85549f1f8d7a
function solving_migration_Kcentric(p, par_def)
    f_tx(t,p,x) = r_Kcentric(x, t, p, par_def)
    t₀ = 0.0
    xspan = (0.0, par_def.col.L)
    prob_tx = ODEProblem(f_tx, t₀, xspan, p)
    solution_tx = solve(prob_tx, alg=par_def.opt.alg, abstol=par_def.opt.abstol, reltol=par_def.opt.reltol)
	tR = solution_tx.u[end]
    return tR
end

# ╔═╡ 8cf3f13b-4c9d-4288-9be8-d53b0dfcb4ad
md"""
### Loss function
"""

# ╔═╡ b2eb2164-1114-4a77-b6a9-37c8945d31ff
function loss_Kcentric(tR, p, col_def, prog_def, sub_0, opt_def)
	# prog_def -> Array of Programs for the different retention times
	# sub_0 -> placeholder Array (1-element) of Substances with 
    tR_p = Array{Any}(undef, length(tR))
    for i=1:length(tR)
		par_def_i = GasChromatographySimulator.Parameters(col_def, prog_def[i], [sub_0], opt_def)
        tR_p[i] = solving_migration_Kcentric(p, par_def_i)
    end
    return sum((tR.-tR_p).^2), tR_p
end

# ╔═╡ 6de1949a-ccf9-4cec-b717-684063787192
measurements[4,:]

# ╔═╡ 6baa0c1c-8b53-4a18-81f2-f2b292f6fcdc
[measurements[4,2], measurements[4,3], measurements[4,4], measurements[4,5]]

# ╔═╡ 1f29fcb5-304b-4e8a-9090-32f877f16dd1
solutes

# ╔═╡ 31229dd0-522d-4b49-8eb9-2a46255435cc
md"""
## Plot of the loss function
"""

# ╔═╡ 199f10e6-2dfd-4fb1-9bf2-8b79f0023091
begin
	Tchar = 320.0:8.0:360.0
	θchar = 20.0:2.0:40.0
	ΔCp = 0.0:20.0:200.0
	retention_time = measurements.tR_Octane[4]
	
	SS = Array{Float64}(undef, length(Tchar), length(θchar), length(ΔCp))
	for k=1:length(ΔCp)
	    for j=1:length(θchar)
	        for i=1:length(Tchar)
	            pp = [Tchar[i], θchar[j], ΔCp[k]]
				sub_0 = sub_def[1]
	            SS[i,j,k] = loss_Kcentric(retention_time, pp, col_def, prog_def, sub_0, opt_def)[1]
	        end
	    end
	end
end

# ╔═╡ 58aa2d2d-cc91-4eea-bc45-3023c61f068a
begin
	plotly()
	plot(plot(Tchar, θchar, SS[:,:,5]', xlabel="Tchar", ylabel="θchar", st=:surface), plot(Tchar, ΔCp, SS[:,5,:]', xlabel="Tchar", ylabel="ΔCp", st=:surface))
end

# ╔═╡ 8eca087c-cb38-4c60-a574-f1161833af53
md"""
### Note
relative flat/shallow minimum
"""

# ╔═╡ 49a8b09f-49e9-49fa-8aef-c9e7bb32f59c
md"""
## Optimization Optim.jl

For now, just use the NelderMead() algorithm. 
"""

# ╔═╡ 78d8789f-952d-4abd-ba73-e540ad164e36
retention_times = [measurements[:,2] measurements[:,3] measurements[:,4] measurements[:,5]]

# ╔═╡ 31f05419-88f6-4dbb-9881-812864b1717f
retention_times_randn = [measurements_randn[:,2] measurements_randn[:,3] measurements_randn[:,4] measurements_randn[:,5]]

# ╔═╡ 4161aef4-0b88-41dd-9f70-925989e2fa8a
f_opt_Kcentric(x) = loss_Kcentric(retention_times[:,1], x, col_def, prog_def, sub_def[1], opt_def)[1]

# ╔═╡ e4c126f3-e705-4f56-9443-026a708ccb90
f_opt_Kcentric_randn(x) = loss_Kcentric(retention_times_randn[:,1], x, col_def, prog_def, sub_def[1], opt_def)[1]

# ╔═╡ 155ff9ce-b6c0-408f-addb-25031832623e
p0_Kcentric = [340.0, 25.0, 50.0]

# ╔═╡ 78a67009-b5bf-4dfd-8397-bf195bc00df0
res_opt_Kcentric_NM = Optim.optimize(f_opt_Kcentric, p0_Kcentric, NelderMead(), Optim.Options(iterations=1000))

# ╔═╡ c8f4ca32-c7de-4087-aa9a-b67d891b220d
res_opt_Kcentric_randn_NM = Optim.optimize(f_opt_Kcentric_randn, p0_Kcentric, NelderMead(), Optim.Options(iterations=1000))

# ╔═╡ a0521a79-cfc6-4c04-80e9-6bf61de778db
Optim.minimizer(res_opt_Kcentric_NM)

# ╔═╡ c70ee0be-5377-47ad-b9b3-82f9d8d7c395
md"""
With the 'exact' retention times and a initial guess near the true value, the optimizde parameters are near the true values.
"""

# ╔═╡ c00b61b7-a5ab-402b-9369-1e6e94b13e68
Optim.minimizer(res_opt_Kcentric_randn_NM)

# ╔═╡ 49a3d205-646a-4796-a07a-c7dbdd0bb720
md"""
With randomly small changes of the retention times, the optimization process is unsucessful and aborts after the maximum of iterations (1000). The parameters Tchar and θchar are near their true value, while the third parameter ΔCp is far off.
"""

# ╔═╡ 36ff90b6-40a1-4dfb-ba9a-e50f6586b474
sub_def[1]

# ╔═╡ 6c241f93-fef6-43d9-aa66-9b385d991ade
begin#=
	Tchar_Optim = Array{Float64}(undef, length(solutes))
	θchar_Optim = Array{Float64}(undef, length(solutes))
	ΔCp_Optim = Array{Float64}(undef, length(solutes))
	res_NM = Array{Any}(undef, length(solutes))
	for i=1:length(solutes)
		f_opt(x) = loss_Kcentric(retention_times[:,i], x, col_def, prog_def, sub_def[i], opt_def)[1]
		p0 = [360.0, 30.0, 100.0]
		res_NM[i] = optimize(f_opt, p0, NelderMead())
		Tchar_Optim[i], θchar_Optim[i], ΔCp_Optim[i] = Optim.minimizer(res_NM[i])
	end
=#end

# ╔═╡ 2fbf3456-85ce-43b0-9d5b-0d0dbb6ec72c
begin#=
	Tchar_Optim_randn = Array{Float64}(undef, length(solutes))
	θchar_Optim_randn = Array{Float64}(undef, length(solutes))
	ΔCp_Optim_randn = Array{Float64}(undef, length(solutes))
	res_randn_NM = Array{Any}(undef, length(solutes))
	for i=1:length(solutes)
		f_opt(x) = loss_Kcentric(retention_times_randn[:,i], x, col_def, prog_def, sub_def[i], opt_def)[1]
		p0 = [360.0, 30.0, 100.0]
		res_randn_NM[i] = optimize(f_opt, p0, NelderMead())
		Tchar_Optim_randn[i], θchar_Optim_randn[i], ΔCp_Optim_randn[i] = Optim.minimizer(res_randn_NM[i])
	end
=#end

# ╔═╡ a704d649-d76e-44f8-81f7-9d36d2cd6a19
sub_def

# ╔═╡ aaf2d59a-71b8-4bd2-96bf-97b70cac44a2
Tchar_Optim

# ╔═╡ b921bebf-883f-4cb8-a814-9fd799ba4531
θchar_Optim

# ╔═╡ 6f0c8b1f-a269-4c7f-a2bf-be4073e91962
ΔCp_Optim

# ╔═╡ 449e3a65-25d7-4274-b067-ec3e4b40d497
measurements

# ╔═╡ 57c08fe2-0fc6-4031-941e-b25dd93687e0
res_NM

# ╔═╡ be32bd9f-13a7-4803-b67f-a484fe41d41c
res_randn_NM

# ╔═╡ b8e40efe-386b-4d08-ae7f-3eb20eeb6865
md"""
### Note
The NelderMead() algorithm with Optim.jl is slow and not fully succesfull (reaches maximum number of itereations). A good estimation of the parameters help?
"""

# ╔═╡ 56a5ac24-b5f6-4136-80c5-6e8928fc1782
md"""
### Simulate with estimated parameters
"""

# ╔═╡ c8269ba6-c06b-4b9d-b422-7db908fa6e40
sub_Octane = GasChromatographySimulator.Substance("Octane", "111-65-9", Tchar_Optim[1], θchar_Optim[1], ΔCp_Optim[1], 0.001, "NelderMead", 0.000107804, 0.0, 0.0)

# ╔═╡ 24cd214c-d984-42b9-8c47-14ffee6ad03e
sub_2Octanone = GasChromatographySimulator.Substance("2-Octanone", "111-13-7", Tchar_Optim[2], θchar_Optim[2], ΔCp_Optim[2], 0.001, "NelderMead", 0.000107102, 0.0, 0.0)

# ╔═╡ d4f9ef63-aca5-4be3-9945-161583947dc6
sub_1Octanol = GasChromatographySimulator.Substance("1-Octanol", "111-87-5", Tchar_Optim[3], θchar_Optim[3], ΔCp_Optim[3], 0.001, "NelderMead", 0.000105558, 0.0, 0.0)

# ╔═╡ 4c62e3af-11de-4d2b-b84f-68cbaa512a25
sub_2Octanol = GasChromatographySimulator.Substance("2-Octanol", "123-96-6", Tchar_Optim[4], θchar_Optim[4], ΔCp_Optim[4], 0.001, "NelderMead", 0.000105558, 0.0, 0.0)

# ╔═╡ 645a76cc-d190-4260-9938-77a2e44a4e79
sub_est = [sub_Octane, sub_2Octanone, sub_1Octanol, sub_2Octanol]

# ╔═╡ c6dcb4b4-479f-4e93-9146-bb420b2207f8
begin
	par_est = Array{GasChromatographySimulator.Parameters}(undef, length(heating_rates))
	for i=1:length(heating_rates)
		par_est[i] = GasChromatographySimulator.Parameters(col_def, prog_def[i], sub_est, opt_def)
	end
end

# ╔═╡ 2211927e-3434-4d7b-ad25-636dfb56c506
begin
	pl_est = Array{DataFrame}(undef, length(heating_rates))
	tR_est = Array{Float64}(undef, length(heating_rates), length(solutes))
	estimations = DataFrame(heating_rate=heating_rates)
	for i=1:length(heating_rates)
		pl_est[i] = GasChromatographySimulator.simulate(par_est[i])[1]
		for j=1:length(solutes)
			jj = findfirst(pl_est[i].Name.==solutes[j])
			tR_est[i,j] = pl_est[i].tR[jj]
		end
	end
	for j=1:length(solutes)
		estimations[!, Symbol(string("tR_", solutes[j]))] = tR_est[:,j]
	end
	estimations
end

# ╔═╡ 8b8c6f52-1ac0-495d-b33b-af75e622f67d
measurements

# ╔═╡ 44a851e4-68da-40f7-9810-0e0fbfb2116c
comparison = DataFrame(heating_rate=measurements.heating_rate, ΔtR_Octane=measurements[:,2]-estimations[:,2],ΔtR_2Octanone=measurements[:,3]-estimations[:,3],ΔtR_1Octanol=measurements[:,4]-estimations[:,4],ΔtR_2Octanol=measurements[:,5]-estimations[:,5])

# ╔═╡ fbb263ea-3018-4c6c-9c49-1f12178e77c5
measurements[:,2]

# ╔═╡ eedc9fdb-6c5a-441a-a437-e08864e2a202
md"""
## Gradient
"""

# ╔═╡ 0dda1a03-4b9b-4e71-bee1-4c9a34a31d4e
g(x) = ForwardDiff.gradient(f_opt_Kcentric, x)

# ╔═╡ cbeed1d3-310f-4717-9015-98ed6a9f092c
g([340.0, 25.0, 50.0])

# ╔═╡ 7727d394-3eed-4eaa-a3bd-8c6b72b1a3c8
sub_def

# ╔═╡ 34c8df74-d931-4d95-891f-8ffa0930abd2
begin
	Tchar0 = sub_def[1].Tchar
	θchar0 = sub_def[1].θchar
	ΔCp0 = sub_def[1].ΔCp
	Tchar_ = 320.0:1.0:360.0
	θchar_ = 20.0:0.1:30.0
	ΔCp_ = 20.0:1.0:100.0
	g_ = Array{Float64}(undef, length(Tchar_))
	g__ = Array{Float64}(undef, length(θchar_))
	g___ = Array{Float64}(undef, length(ΔCp_))
	for i=1:length(Tchar_)
		g_[i] = g([Tchar_[i], θchar0, ΔCp0])[1]
	end
	for i=1:length(θchar_)
		g__[i] = g([Tchar0, θchar_[i], ΔCp0])[2]
	end
	for i=1:length(ΔCp_)
		g___[i] = g([Tchar0, θchar0, ΔCp_[i]])[3]
	end
end

# ╔═╡ 4b3bf6b1-2046-406f-953b-47fbe1505dd0
plot(Tchar_, g_)

# ╔═╡ 8bbff1cd-7e66-4510-9af1-4271aeb52be2
plot(θchar_, g__)

# ╔═╡ 303f07c5-ee4b-470a-9e16-51786953dde0
plot(ΔCp_, g___)

# ╔═╡ e6198bf4-13d5-48d6-94c8-2b59d9e0932f
md"""
### Use IntervalRootFinding.jl on the gradient
"""

# ╔═╡ cd7936bd-cdd6-4a67-ae46-01da72bd5c80
gg((x,y,z)) = SVector{3}(g([x,y,z]))

# ╔═╡ 648ac735-fa33-4bd2-af09-b3dca760f4da
gg((340.0, 25.0, 60.0))

# ╔═╡ 5438631f-f714-4121-ae7e-334bcebb52b9
md"""
Unknown Error
"""

# ╔═╡ 06457aae-9664-4e97-ba49-c19af869ed44
roots(gg, (320.0..360.0) × (25.0..30.0) × (50.0..100.0), Krawczyk)

# ╔═╡ 77261d7a-4d97-4a33-89ca-593a805c72dd
function G( (x1, x2, x3) )
           return SVector(x1^2 + x2^2 + x3^2 - 1,
                          x1^2 + x3^2 - 0.25,
                          x1^2 + x2^2 - 4x3
                         )
       end

# ╔═╡ edb976ae-10b6-4d7c-824e-c2e53a438697
X = -5..5

# ╔═╡ 5fe994c2-0415-41c6-b496-c6d2605b141a
G((1, 3, -1))

# ╔═╡ d1096159-a9f0-44ab-82a1-baffe0f77b04
rts = roots(G, X × X × X)

# ╔═╡ 5733f75a-f2da-40b9-a0e9-4fa95e0ea5bb
#res_opt_Kcentric_BFGS = optimize(f_opt_Kcentric, p0_Kcentric, BFGS())

# ╔═╡ 3164d218-3d1a-48dd-9edf-e6b9d644b4b3
md"""
## 👎 Optimization LeastSquaresOptim.jl
"""

# ╔═╡ 21db8933-4ef8-4f47-a9c8-a725960b7e50
p0_Kcentric

# ╔═╡ 8ea1e90a-1c1c-4c1d-b464-df389e0da82c
md"""
All working combinations have parameters Tchar and θchar near their true value, BUT ΔCp significantly lower.
"""

# ╔═╡ 6155b441-2582-4d95-935c-0ff277b05e28
md"""
### LevenbergMarquardt()
"""

# ╔═╡ 631478a6-9956-4627-a83e-4caa0d619bda
LeastSquaresOptim.optimize(f_opt_Kcentric, p0_Kcentric, LevenbergMarquardt())

# ╔═╡ 52991ab7-80c6-4b52-9f49-85af71f48667
LeastSquaresOptim.optimize(f_opt_Kcentric, p0_Kcentric, LevenbergMarquardt(LeastSquaresOptim.Cholesky()))

# ╔═╡ 67d9ed51-b60c-4b77-a61e-37e7b649e06b
LeastSquaresOptim.optimize(f_opt_Kcentric, p0_Kcentric, LevenbergMarquardt(LeastSquaresOptim.LSMR()))

# ╔═╡ 0e002313-aadc-4308-ba75-4e55160eb538
# try different initial guesses

# ╔═╡ e25a94dd-227b-46db-b8fa-b1596704c0d5
sub_def[1]

# ╔═╡ c74db1b3-6e43-4cf8-b799-f134eeb9d58d
LeastSquaresOptim.optimize(f_opt_Kcentric, [340.5, 26.8, 66.4], LevenbergMarquardt())

# ╔═╡ a95cfdd2-bb43-4775-a166-26999d3380d9
#LeastSquaresOptim.optimize(f_opt_Kcentric, [300.0, 20.0, 40.0], LevenbergMarquardt())

# ╔═╡ 011c0ef3-0321-43e4-aabd-5dc79021f10b
md"""
### Dogleg()
"""

# ╔═╡ e943d066-afb0-4ba7-9dc7-8f6e1b48879f
LeastSquaresOptim.optimize(f_opt_Kcentric_randn, p0_Kcentric, Dogleg())

# ╔═╡ a238f44c-f923-47b6-9ed4-e52c23b591fe
LeastSquaresOptim.optimize(f_opt_Kcentric_randn, p0_Kcentric, Dogleg(LeastSquaresOptim.QR()))

# ╔═╡ 33665a61-1169-4536-b6e3-b102ca3dd0ea
LeastSquaresOptim.optimize(f_opt_Kcentric_randn, p0_Kcentric, Dogleg(LeastSquaresOptim.Cholesky()))

# ╔═╡ e0b3f4e7-1457-4518-928a-22c5fda17384
LeastSquaresOptim.optimize(f_opt_Kcentric_randn, p0_Kcentric, Dogleg(LeastSquaresOptim.LSMR()))

# ╔═╡ e0233619-1766-4e4c-b23f-abdae898ab1c
sub_def[1]

# ╔═╡ 7fe26669-bcdc-4242-a1da-f4046a983e55
md"""
## BlackBoxOptim.jl

Slow, but good results
"""

# ╔═╡ 7f1341eb-e348-460a-b945-25ffa691ce06
res_bb = bboptimize(f_opt_Kcentric; SearchRange=[(300.0, 360.0),(20.0, 40.0), (50.0, 100.0)])

# ╔═╡ 05dd1d4a-1b8e-4a43-afe2-5dd358f3e6d0
# stopped after 10000 iteration but with a good fit

# ╔═╡ a83db724-2c91-4865-8aa0-58a5c679d769
sub_def[1]

# ╔═╡ 720eecfa-03dd-4412-a974-7208e83c5841
res_randn_bb = bboptimize(f_opt_Kcentric_randn; SearchRange=[(300.0, 360.0),(20.0, 40.0), (50.0, 100.0)])

# ╔═╡ e9d6b0e6-939b-497e-8a16-58253a73d4ba
comp = compare_optimizers(f_opt_Kcentric; SearchRange=[(300.0, 360.0),(20.0, 40.0), (50.0, 100.0)], MaxTime = 10.0)

# ╔═╡ a21b634c-2c81-442f-bfb0-15f9e101605d
comp

# ╔═╡ f26ac00e-a12e-4399-a52b-528202e859c6
#initial guess:
begin
	Tchar_guess = 340.0
	θchar_guess = 22.0*(Tchar_guess/273.15)^0.7
	ΔCp_guess = 100.0
	guess = [Tchar_guess, θchar_guess, ΔCp_guess]
	res_bb_guess = bboptimize(f_opt_Kcentric, guess; SearchRange=[(300.0, 360.0),(20.0, 40.0), (50.0, 100.0)])
end

# ╔═╡ 28ce4c23-7956-45d8-911d-e5ea64e23c3b
# good fits from the compare_optimizers:
# :de_rand_1_bin_radiuslimited
# :resampling_memetic_search

# ╔═╡ b0ec794c-3b6b-4395-8e77-f75c2d78b2af
res_bb_guess_1 = bboptimize(f_opt_Kcentric, guess; SearchRange=[(300.0, 360.0),(20.0, 40.0), (50.0, 100.0)], Method=:de_rand_1_bin_radiuslimited)

# ╔═╡ 045c3a97-c2fc-4637-806f-0a621c33dc43
res_bb_guess_small_range = bboptimize(f_opt_Kcentric, [guess[1], guess[2], 65.0]; SearchRange=[(335.0, 345.0),(22.0, 28.0), (60.0, 70.0)])

# ╔═╡ 6b0b29d4-276a-4888-8d17-a4468b98af95
md"""
# End
"""

# ╔═╡ Cell order:
# ╠═bdaf89aa-42cd-4ed4-b009-a7820e233cc9
# ╠═2c5e497e-09f2-4ffe-8df5-c455fd275771
# ╠═f3ab4008-8851-4ba2-87bc-eb54ebbc0d7d
# ╠═02ec232f-83f8-4b60-afed-904cdc08e106
# ╠═33d48cb9-a06c-4669-9df3-a93de198e8a9
# ╠═fdecf17b-0b87-4df0-9bb5-a946505ce290
# ╠═6039f94d-baf2-4757-9779-29791b1addc9
# ╠═d9473aa5-40cd-47aa-ad5a-76c9ebfc1ade
# ╠═4289f252-d7ba-4a21-aea3-1555e3d22d22
# ╠═dab67cfa-dd79-44d4-bd78-327dafb550ed
# ╠═8bc58767-2d14-4317-b94b-8c47e2564f30
# ╠═c2a06f83-0c27-4841-ba3e-fd8eb34565a2
# ╠═480a7a33-d13c-40ce-83ee-28ff647e7a70
# ╠═e8368d40-faae-4a4c-b88a-985da383a4c0
# ╠═5f9870a7-28c8-49a8-9ba6-77b94df6c908
# ╠═aac1d065-18a5-4109-8dc2-510b3a6cab0c
# ╠═d95f2b44-f7b7-4d6c-a9c8-7b2d5b462b37
# ╠═c22086bf-17c5-4608-8e44-26a30f0689d2
# ╠═1cd045bd-4d72-4b4d-94be-ceb452f7e989
# ╠═c5b85477-ed0a-42b7-bd4e-aa1282871b19
# ╟─4fcc7a44-7899-45ac-a4d6-b1e53bec411d
# ╠═708d03c2-546e-4f89-be85-f9d57e656a01
# ╠═94424a15-c0fb-4f0c-be1c-6146585f1e25
# ╠═4919d447-709e-4261-acb2-e713937f86fc
# ╠═00603cca-e1bd-4db7-af0f-a7661edbc109
# ╠═34dd3b26-4ef8-47da-837a-a6daee7356d8
# ╠═9b9f492b-d4aa-4ac1-9906-51b3f29ce17d
# ╠═9ccf93d6-0426-4d33-86f2-21537811bfcb
# ╠═98cb52b0-65a4-4fde-836b-e8837363c00e
# ╠═ed18e731-67ba-433e-a802-a258f097bca4
# ╠═4d144797-380b-45e0-b412-63b08bd21762
# ╠═f48310e7-8ffc-43cd-be23-169c418ba80d
# ╠═906da7f5-ba06-4ecd-9850-084a938cafc4
# ╠═915a8459-2265-4a7d-9113-de562dfc02a5
# ╠═322d274a-bfc7-46af-999a-6d3076f905c3
# ╠═2446c0d3-8105-4644-a67d-ec39ec4dc1fc
# ╠═5b136b18-5a0f-4cff-8a1c-585b390bab0c
# ╠═f1c1d843-ef16-4212-81c1-56ede0301c54
# ╠═5ce58c17-ec49-4fa1-abed-961cb1380d6a
# ╠═599cf329-ed37-453c-a190-b909847c1ec9
# ╠═7011dff6-2e07-4d03-b21f-94fd4359adc6
# ╠═9fbd5c15-1bb5-4344-b6a4-85549f1f8d7a
# ╠═8cf3f13b-4c9d-4288-9be8-d53b0dfcb4ad
# ╠═b2eb2164-1114-4a77-b6a9-37c8945d31ff
# ╠═6de1949a-ccf9-4cec-b717-684063787192
# ╠═6baa0c1c-8b53-4a18-81f2-f2b292f6fcdc
# ╠═1f29fcb5-304b-4e8a-9090-32f877f16dd1
# ╠═31229dd0-522d-4b49-8eb9-2a46255435cc
# ╠═199f10e6-2dfd-4fb1-9bf2-8b79f0023091
# ╠═58aa2d2d-cc91-4eea-bc45-3023c61f068a
# ╠═8eca087c-cb38-4c60-a574-f1161833af53
# ╠═49a8b09f-49e9-49fa-8aef-c9e7bb32f59c
# ╠═78d8789f-952d-4abd-ba73-e540ad164e36
# ╠═31f05419-88f6-4dbb-9881-812864b1717f
# ╠═4161aef4-0b88-41dd-9f70-925989e2fa8a
# ╠═e4c126f3-e705-4f56-9443-026a708ccb90
# ╠═155ff9ce-b6c0-408f-addb-25031832623e
# ╠═78a67009-b5bf-4dfd-8397-bf195bc00df0
# ╠═c8f4ca32-c7de-4087-aa9a-b67d891b220d
# ╠═a0521a79-cfc6-4c04-80e9-6bf61de778db
# ╟─c70ee0be-5377-47ad-b9b3-82f9d8d7c395
# ╠═c00b61b7-a5ab-402b-9369-1e6e94b13e68
# ╟─49a3d205-646a-4796-a07a-c7dbdd0bb720
# ╠═36ff90b6-40a1-4dfb-ba9a-e50f6586b474
# ╠═6c241f93-fef6-43d9-aa66-9b385d991ade
# ╠═2fbf3456-85ce-43b0-9d5b-0d0dbb6ec72c
# ╠═a704d649-d76e-44f8-81f7-9d36d2cd6a19
# ╠═aaf2d59a-71b8-4bd2-96bf-97b70cac44a2
# ╠═b921bebf-883f-4cb8-a814-9fd799ba4531
# ╠═6f0c8b1f-a269-4c7f-a2bf-be4073e91962
# ╠═449e3a65-25d7-4274-b067-ec3e4b40d497
# ╠═57c08fe2-0fc6-4031-941e-b25dd93687e0
# ╠═be32bd9f-13a7-4803-b67f-a484fe41d41c
# ╠═b8e40efe-386b-4d08-ae7f-3eb20eeb6865
# ╠═56a5ac24-b5f6-4136-80c5-6e8928fc1782
# ╠═c8269ba6-c06b-4b9d-b422-7db908fa6e40
# ╠═24cd214c-d984-42b9-8c47-14ffee6ad03e
# ╠═d4f9ef63-aca5-4be3-9945-161583947dc6
# ╠═4c62e3af-11de-4d2b-b84f-68cbaa512a25
# ╠═645a76cc-d190-4260-9938-77a2e44a4e79
# ╠═c6dcb4b4-479f-4e93-9146-bb420b2207f8
# ╠═2211927e-3434-4d7b-ad25-636dfb56c506
# ╠═8b8c6f52-1ac0-495d-b33b-af75e622f67d
# ╠═44a851e4-68da-40f7-9810-0e0fbfb2116c
# ╠═fbb263ea-3018-4c6c-9c49-1f12178e77c5
# ╠═eedc9fdb-6c5a-441a-a437-e08864e2a202
# ╠═0dda1a03-4b9b-4e71-bee1-4c9a34a31d4e
# ╠═cbeed1d3-310f-4717-9015-98ed6a9f092c
# ╠═7727d394-3eed-4eaa-a3bd-8c6b72b1a3c8
# ╠═34c8df74-d931-4d95-891f-8ffa0930abd2
# ╠═4b3bf6b1-2046-406f-953b-47fbe1505dd0
# ╠═8bbff1cd-7e66-4510-9af1-4271aeb52be2
# ╠═303f07c5-ee4b-470a-9e16-51786953dde0
# ╠═e6198bf4-13d5-48d6-94c8-2b59d9e0932f
# ╠═cd7936bd-cdd6-4a67-ae46-01da72bd5c80
# ╠═648ac735-fa33-4bd2-af09-b3dca760f4da
# ╠═5438631f-f714-4121-ae7e-334bcebb52b9
# ╠═06457aae-9664-4e97-ba49-c19af869ed44
# ╠═77261d7a-4d97-4a33-89ca-593a805c72dd
# ╠═edb976ae-10b6-4d7c-824e-c2e53a438697
# ╠═5fe994c2-0415-41c6-b496-c6d2605b141a
# ╠═d1096159-a9f0-44ab-82a1-baffe0f77b04
# ╠═5733f75a-f2da-40b9-a0e9-4fa95e0ea5bb
# ╠═3164d218-3d1a-48dd-9edf-e6b9d644b4b3
# ╠═21db8933-4ef8-4f47-a9c8-a725960b7e50
# ╠═8ea1e90a-1c1c-4c1d-b464-df389e0da82c
# ╠═6155b441-2582-4d95-935c-0ff277b05e28
# ╠═631478a6-9956-4627-a83e-4caa0d619bda
# ╠═52991ab7-80c6-4b52-9f49-85af71f48667
# ╠═67d9ed51-b60c-4b77-a61e-37e7b649e06b
# ╠═0e002313-aadc-4308-ba75-4e55160eb538
# ╠═e25a94dd-227b-46db-b8fa-b1596704c0d5
# ╠═c74db1b3-6e43-4cf8-b799-f134eeb9d58d
# ╠═a95cfdd2-bb43-4775-a166-26999d3380d9
# ╠═011c0ef3-0321-43e4-aabd-5dc79021f10b
# ╠═e943d066-afb0-4ba7-9dc7-8f6e1b48879f
# ╠═a238f44c-f923-47b6-9ed4-e52c23b591fe
# ╠═33665a61-1169-4536-b6e3-b102ca3dd0ea
# ╠═e0b3f4e7-1457-4518-928a-22c5fda17384
# ╠═e0233619-1766-4e4c-b23f-abdae898ab1c
# ╠═7fe26669-bcdc-4242-a1da-f4046a983e55
# ╠═7f1341eb-e348-460a-b945-25ffa691ce06
# ╠═05dd1d4a-1b8e-4a43-afe2-5dd358f3e6d0
# ╠═a83db724-2c91-4865-8aa0-58a5c679d769
# ╠═720eecfa-03dd-4412-a974-7208e83c5841
# ╠═e9d6b0e6-939b-497e-8a16-58253a73d4ba
# ╠═a21b634c-2c81-442f-bfb0-15f9e101605d
# ╠═f26ac00e-a12e-4399-a52b-528202e859c6
# ╠═28ce4c23-7956-45d8-911d-e5ea64e23c3b
# ╠═b0ec794c-3b6b-4395-8e77-f75c2d78b2af
# ╠═045c3a97-c2fc-4637-806f-0a621c33dc43
# ╠═6b0b29d4-276a-4888-8d17-a4468b98af95
