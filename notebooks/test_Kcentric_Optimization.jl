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
	using Plots
	using Optimization
	using OptimizationOptimJL
	using OptimizationBBO
	using Zygote
	using DiffEqSensitivity
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

# ╔═╡ 5ce58c17-ec49-4fa1-abed-961cb1380d6a
pl_meas

# ╔═╡ 78d8789f-952d-4abd-ba73-e540ad164e36
retention_times = [measurements[:,2] measurements[:,3] measurements[:,4] measurements[:,5]]

# ╔═╡ 31f05419-88f6-4dbb-9881-812864b1717f
retention_times_randn = [measurements_randn[:,2] measurements_randn[:,3] measurements_randn[:,4] measurements_randn[:,5]]

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

# ╔═╡ 1574adf8-9abe-49e3-8bfd-9970009a9373
opt_def_ = GasChromatographySimulator.Options(ng=true) # setting ng=true

# ╔═╡ 4f1b60f1-0a08-4983-a051-b278afbea4bb
# function with modified equation for retention factor k(x,t) -> no correction for dimless film thickness φ changed in realtion to the dimless film thickness used for the estimation of parameters φ₀ -> here both are the same
function tR_calc(Tchar, θchar, ΔCp, L, d, prog, opt, gas)
	# df has no influence on the result (is hidden in Tchar, θchar, ΔCp)
	R = 8.3145
	k(x,t,Tchar,θchar,ΔCp) = exp((ΔCp/R + Tchar/θchar)*(Tchar/prog.T_itp(x,t)-1) + ΔCp/R*log(prog.T_itp(x,t)/Tchar))
	rM(x,t,L,d) = GasChromatographySimulator.mobile_phase_residency(x,t, prog.T_itp, prog.Fpin_itp, prog.pout_itp, L, d, gas; ng=opt.ng, vis=opt.vis, control=opt.control)
	r(t,p,x) = (1+k(x,t,p[1],p[2],p[3]))*rM(x,t,p[4],p[5])
	t₀ = 0.0
	xspan = (0.0, L)
	p = [Tchar, θchar, ΔCp, L, d]
	prob = ODEProblem(r, t₀, xspan, p)
	solution = solve(prob, alg=opt.alg, abstol=opt.abstol, reltol=opt.reltol)
	tR = solution.u[end]
	return tR
end

# ╔═╡ 15f73a94-8d6f-4b1d-8314-076f8022f302
Tchar0 = [340.0, 380.0, 380.0, 390.0]

# ╔═╡ b72343d0-e831-4699-9dcb-5251f4cf5738
θchar0 = 22.0.*(Tchar0./273.15).^0.7

# ╔═╡ e20aeeed-5815-484c-acd7-ea57598b3642
ΔCp0 = [70.0, 80.0, 90.0, 90.0]

# ╔═╡ 03f40097-f82b-4533-a3b4-9fbe3c3d67b3
L0 = 29.8

# ╔═╡ 79d87deb-bf52-4d39-aae1-b16040d3ac62
d0 = 0.26e-3

# ╔═╡ 12cd49e8-9ed5-47e8-8415-5af39c49faf1
p0 = [Tchar0, θchar0, ΔCp0, L0, d0]

# ╔═╡ 829a7d93-e963-4084-b691-d9c12b41d0e2
p0[1]

# ╔═╡ ad81e670-eec3-4610-92b7-99da8a9eb658
p0__ = [Tchar0; θchar0; ΔCp0]

# ╔═╡ 3176cd9f-1c63-4b69-ab74-85a01a24f54a
size(retention_times)

# ╔═╡ 8cf3f13b-4c9d-4288-9be8-d53b0dfcb4ad
md"""
### Loss function
"""

# ╔═╡ d002c5fb-43b1-4873-8245-17fc94cde9c3
function loss(tR, p, prog, opt)
	# loss function as sum over all programs and solutes
	gas = "He"
	ns = size(tR)[2]
	Tchar0 = p[1:ns] # Array length = number solutes
	θchar0 = p[ns+1:2*ns] # Array length = number solutes
	ΔCp0 = p[2*ns+1:3*ns] # Array length = number solutes
	L0 = p[end-1] # Float64
	d0 = p[end] # Float64
	tRcalc = Array{Float64}(undef, size(tR)[1], size(tR)[2])
	for j=1:size(tR)[2]
		for i=1:size(tR)[1]
			tRcalc[i,j] = tR_calc(Tchar0[j], θchar0[j], ΔCp0[j], L0, d0, prog[i], opt, gas)
		end
	end
	return sum((tR.-tRcalc).^2), tRcalc
end	

# ╔═╡ 42cbd777-f727-492d-b04f-49296f6be31f
function loss(tR, p, L, d, prog, opt)
	# loss function as sum over all programs and solutes, fixed L and d
	gas = "He"
	ns = size(tR)[2]
	Tchar0 = p[1:ns] # Array length = number solutes
	θchar0 = p[ns+1:2*ns] # Array length = number solutes
	ΔCp0 = p[2*ns+1:3*ns] # Array length = number solutes
	#L0 = p[end-1] # Float64
	#d0 = p[end] # Float64
	tRcalc = Array{Float64}(undef, size(tR)[1], size(tR)[2])
	for j=1:size(tR)[2]
		for i=1:size(tR)[1]
			tRcalc[i,j] = tR_calc(Tchar0[j], θchar0[j], ΔCp0[j], L, d, prog[i], opt, gas)
		end
	end
	return sum((tR.-tRcalc).^2), tRcalc
end	

# ╔═╡ 29ab9ff2-7472-4a95-80c8-6d64e4803a42
function loss(p_opt, p)
	# loss function as sum over all programs and solutes, fixed L and d
	gas = "He"
	tR = p[1]
	L = p[2]
	d = p[3]
	prog = p[4]
	opt = p[5]
	
	ns = size(tR)[2]
	Tchar0 = p_opt[1:ns] # Array length = number solutes
	θchar0 = p_opt[ns+1:2*ns] # Array length = number solutes
	ΔCp0 = p_opt[2*ns+1:3*ns] # Array length = number solutes
	#L0 = p[end-1] # Float64
	#d0 = p[end] # Float64
	tRcalc = Array{Float64}(undef, size(tR)[1], size(tR)[2])
	for j=1:size(tR)[2]
		for i=1:size(tR)[1]
			tRcalc[i,j] = tR_calc(Tchar0[j], θchar0[j], ΔCp0[j], L, d, prog[i], opt, gas)
		end
	end
	return sum((tR.-tRcalc).^2), tRcalc
end	

# ╔═╡ ef392841-68b5-4407-8920-954cd4d78c4c
function loss_(tR, p, prog, opt)
	# loss function as sum over all programs (for one solute)
	gas = "He"
	Tchar0 = p[1] # Float64
	θchar0 = p[2] # Float64
	ΔCp0 = p[3] # Float64
	L0 = p[4] # Float64
	d0 = p[5] # Float64
	tRcalc = Array{Float64}(undef, length(tR))
	for i=1:length(tR)
		tRcalc[i] = tR_calc(Tchar0, θchar0, ΔCp0, L0, d0, prog[i], opt, gas)
	end
	return sum((tR.-tRcalc).^2), tRcalc
end	

# ╔═╡ 71d799f2-b7e0-4145-9772-2fb13b0400a7
loss(retention_times, p0__, prog_def, opt_def_)

# ╔═╡ 672bd169-d9e7-4825-aa83-733d71e347d5
p0_ = [Tchar0[1], θchar0[1], ΔCp0[1], L0, d0]

# ╔═╡ 6510d72a-01d7-46de-9385-e76b6167a8b2
loss_(retention_times[:,1], p0_, prog_def, opt_def_)

# ╔═╡ 31229dd0-522d-4b49-8eb9-2a46255435cc
md"""
## Plot of the loss function
"""

# ╔═╡ 199f10e6-2dfd-4fb1-9bf2-8b79f0023091
begin
	Tchar = 335.0:1.0:345.0
	θchar = 20.0:1.0:30.0
	ΔCp = 0.0:20.0:200.0
	retention_time = measurements.tR_Octane
	
	SS = Array{Float64}(undef, length(Tchar), length(θchar), length(ΔCp))
	for k=1:length(ΔCp)
	    for j=1:length(θchar)
	        for i=1:length(Tchar)
	            pp = [Tchar[i], θchar[j], ΔCp[k], 30.0, 0.25e-3]
	            SS[i,j,k] = loss_(retention_time, pp, prog_def, opt_def_)[1]
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
## Compile optimization function
"""

# ╔═╡ 4161aef4-0b88-41dd-9f70-925989e2fa8a
f_opt(x,p) = loss(retention_times, x, p[1], p[2], prog_def, opt_def_)[1]

# ╔═╡ e4c126f3-e705-4f56-9443-026a708ccb90
f_opt_randn(x,p) = loss(retention_times_randn, x, p[1], p[2], prog_def, opt_def_)[1]

# ╔═╡ eedc9fdb-6c5a-441a-a437-e08864e2a202
md"""
## Optimization.jl - fixed L, d
"""

# ╔═╡ f948ecbc-f49f-40c7-b697-4a38cb95d458
org = [sub_def[1].Tchar, sub_def[2].Tchar, sub_def[4].Tchar, sub_def[3].Tchar, sub_def[1].θchar, sub_def[2].θchar, sub_def[4].θchar, sub_def[3].θchar, sub_def[1].ΔCp, sub_def[2].ΔCp, sub_def[4].ΔCp, sub_def[3].ΔCp]

# ╔═╡ c7678500-8e16-46de-9583-4d36a0262304
md"""
### OptimizationOptimJL - NelderMead
"""

# ╔═╡ aff1eaa8-c680-407b-8a0f-dadae568db8b
function loss__(p_opt, p)
	# loss function as sum over all programs and solutes, fixed L and d
	gas = "He"
	tR = p[1]
	L = p[2]
	d = p[3]
	prog = p[4]
	opt = p[5]

	if length(size(tR)) == 1
		ns = 1
	else
		ns = size(tR)[2]
	end
	Tchar0 = p_opt[1:ns] # Array length = number solutes
	θchar0 = p_opt[ns+1:2*ns] # Array length = number solutes
	ΔCp0 = p_opt[2*ns+1:3*ns] # Array length = number solutes
	#L0 = p[end-1] # Float64
	#d0 = p[end] # Float64
	tRcalc = Array{Any}(undef, size(tR)[1], ns)
	#buf = Zygote.Buffer(tR, size(tR)[1], ns)
	for j=1:ns
		for i=1:size(tR)[1]
			tRcalc[i,j] = tR_calc(Tchar0[j], θchar0[j], ΔCp0[j], L, d, prog[i], opt, gas)
			#buf[i,j] = tR_calc(Tchar0[j], θchar0[j], ΔCp0[j], L, d, prog[i], opt, gas)
		end
	end
	#tRcalc = copy(buf)
	return sum((tR.-tRcalc).^2), tRcalc
end	

# ╔═╡ 892cf99b-7672-4046-afdd-179673714235
p0__

# ╔═╡ 3784ab35-e63f-40a5-95ad-ff4183519ac6
pp__ = [retention_times, 30.0, 0.25e-3, prog_def, opt_def_]

# ╔═╡ 8aae0422-b0ec-4b61-88f2-1e7b7183a149
loss__(p0__, pp__)

# ╔═╡ 7d0e83ca-4617-4d86-9266-507e0aa31d80
opt_loss__(x, p) = loss__(x, p)[1]

# ╔═╡ 8044d719-1fb1-4a50-87fc-45c5cf8de112
#prob = OptimizationProblem(f_opt, p0__[1:end-2], [30.0, 0.25e-3])
prob = OptimizationProblem(opt_loss__, p0__, pp__)

# ╔═╡ 884f9da5-0f32-423c-a71f-07abcd22f873
sol = solve(prob, NelderMead())

# ╔═╡ 6dab0dcf-c402-4bc8-a6f4-7c7cf305967f
(sol.-org)./org.*100.0

# ╔═╡ b380629a-880d-457d-bb2d-2cc6448ac4f5
md"""
- Good results for ``T_{char}``, <<1% difference
- Ok results for ``\theta_{char}``, <2.5% difference
- not so good results for ``ΔC_p``, up to 11.5% difference
"""

# ╔═╡ 99d7ef10-1e01-4b0d-8bc3-98a761646149
md"""
#### Optimization for only one solute
"""

# ╔═╡ f68804fa-6f32-48cd-9df9-007302928a33
length(size(retention_times[:,1]))

# ╔═╡ 07d0008d-06ee-42ba-8530-9440028d1746
opt_loss__(p0__[1:4:end], [retention_times[:,1], 30.0, 0.25e-3, prog_def, opt_def_])

# ╔═╡ 0af253f5-26ca-4423-84be-ab1f6165f1cf
prob_sol1 = OptimizationProblem(opt_loss__, p0__[1:4:end], [retention_times[:,1], 30.0, 0.25e-3, prog_def, opt_def_])

# ╔═╡ 51f050ce-840e-4a06-af46-bfeb452dc30d
prob_sol2 = OptimizationProblem(opt_loss__, p0__[2:4:end], [retention_times[:,2], 30.0, 0.25e-3, prog_def, opt_def_])

# ╔═╡ abd86076-f11e-4171-9176-21c3779b3bc4
prob_sol3 = OptimizationProblem(opt_loss__, p0__[3:4:end], [retention_times[:,3], 30.0, 0.25e-3, prog_def, opt_def_])

# ╔═╡ 2b509fc6-6ee9-44ca-bb95-df05246f2ca8
prob_sol4 = OptimizationProblem(opt_loss__, p0__[4:4:end], [retention_times[:,4], 30.0, 0.25e-3, prog_def, opt_def_])

# ╔═╡ c685bdc4-a9d0-4cbc-8f33-c7d560c10b5c
sol_sol1 = solve(prob_sol1, NelderMead())

# ╔═╡ 10c316bc-6389-44dd-8a92-b772ffed9c84
sol_sol2 = solve(prob_sol2, NelderMead())

# ╔═╡ 394deaba-a74e-45d0-9f02-c9ec7799e372
sol_sol3 = solve(prob_sol3, NelderMead())

# ╔═╡ 4a909f87-f16b-4b88-bdc7-578e5849f01c
sol_sol4 = solve(prob_sol4, NelderMead())

# ╔═╡ 26bc4084-4a58-4bbe-b0c2-4b087789a5c5
sol_single = [sol_sol1[1], sol_sol2[1], sol_sol3[1], sol_sol4[1], sol_sol1[2], sol_sol2[2], sol_sol3[2], sol_sol4[2], sol_sol1[3], sol_sol2[3], sol_sol3[3], sol_sol4[3]]

# ╔═╡ dcc3a436-d079-47c1-808f-88179fa3d421
org

# ╔═╡ 7481f4bf-6622-495d-84de-d50d02b931f0
(sol_single.-org)./org.*100.0

# ╔═╡ fb68eaa0-ded3-40ac-90ab-ff0d21ed347e
md"""
Different results, if the parameters are optimized for only on solute at a time. Mostly better results:
- Good results for ``T_{char}``, <<0.1% difference
- Ok results for ``\theta_{char}``, <1.6% difference
- not so good results for ``ΔC_p``, up to 14% difference (for one solute worse results)
"""

# ╔═╡ 6d2408b0-2edf-48f5-aa13-83ff8d0c9955
md"""
#### Random modified retention times
"""

# ╔═╡ 46751f67-d7ba-4e8b-8a82-c47becd9fcfd
pp__randn = [retention_times_randn, 30.0, 0.25e-3, prog_def, opt_def_]

# ╔═╡ 9f6122b4-ba52-4029-b86c-84320d7b6969
prob_randn = OptimizationProblem(opt_loss__, p0__, pp__randn)

# ╔═╡ f74ca574-f9d5-4b0f-b190-4a5a27ad78b1
sol_randn = solve(prob_randn, NelderMead())

# ╔═╡ b7cae6f2-ffd7-4760-b3e6-1b869581d3b9
(sol_randn.-org)./org.*100.0

# ╔═╡ f527f9c9-7bab-4fec-b0d0-2e24842a83af
md"""
Similar results for the random modified retention times
- Good results for ``T_{char}``, <<1% difference
- Ok results for ``\theta_{char}``, <4% difference
- not so good results for ``ΔC_p``, up to 14% difference
"""

# ╔═╡ be703ad1-a05a-459f-a437-f02f97cef40b
md"""
### OptimizationBBO
"""

# ╔═╡ 50a81a4a-a51f-4fb9-a104-07416fc479ae
lb = p0__*0.8

# ╔═╡ 3d0f67c6-881e-4a77-b15f-be180f2bb567
ub = p0__*1.2

# ╔═╡ 4e5b0aff-7c6f-46f9-8575-1d4a302973ef
prob_bbo = OptimizationProblem(opt_loss__, p0__, lb=lb, ub=ub, pp__)

# ╔═╡ ab233b0f-62bd-45a0-aaef-d29bcdc7182e
sol_bbo = solve(prob_bbo, BBO_adaptive_de_rand_1_bin_radiuslimited())

# ╔═╡ 7d5bbbf8-6736-4e87-bae1-56e6fe4c4b22
(sol_bbo.-org)./org.*100.0

# ╔═╡ 5a39f6a5-6c1a-4dbe-a6a7-15cbbca37ea8
md"""
- Good results for ``T_{char}``, <<1% difference
- not so good  results for ``\theta_{char}``, mostly <2.0%
- not so good results for ``ΔC_p``, up to 10% difference
"""

# ╔═╡ e4f421e0-516c-493e-ba62-36b93aaee21f
md"""
#### Random modified retention times
"""

# ╔═╡ 6b255417-2e5e-4fe6-9acf-feda9c4713d0
prob_bbo_randn = OptimizationProblem(opt_loss__, p0__, lb=lb, ub=ub, pp__randn)

# ╔═╡ 1c3ec334-4200-473e-a50e-8891f65f2110
sol_bbo_randn = solve(prob_bbo_randn, BBO_adaptive_de_rand_1_bin_radiuslimited());

# ╔═╡ d7682fe5-6297-4220-965a-4a0a833e36b3
(sol_bbo_randn.-org)./org.*100.0

# ╔═╡ 9e923fe7-7059-4c4e-972c-36618ccc5149
md"""
Similar results for the random modified retention times
- Good results for ``T_{char}``, <<1% difference
- not so good  results for ``\theta_{char}``, mostly <2.0% difference, but one with 4% difference
- not so good results for ``ΔC_p``, up to 20% difference
"""

# ╔═╡ d9b5253c-96ac-4d45-a43e-8aabd3f54b49
md"""
### OptimizationOptimJL - BFGS
"""

# ╔═╡ 4fa6d902-8b70-4688-a3c2-f48e42919a66
pp = [retention_times, 30.0, 0.25e-3, prog_def, opt_def_]

# ╔═╡ cd7a1778-0888-432a-b7ba-a5c0e0a94b32
opt_f_fwd = OptimizationFunction(opt_loss__, Optimization.AutoForwardDiff())

# ╔═╡ a2d60fc8-79c2-4c95-90b5-b83ebbf6357e
prob_fwd = OptimizationProblem(opt_f_fwd, p0__, pp, lb=lb, ub=ub)

# ╔═╡ e02f8b49-c372-40fe-80e6-af47f4a75d56
sol_fwd = solve(prob_fwd, BFGS())

# ╔═╡ c362bc14-4f0c-4111-b318-e640c96fb157
(sol_fwd.-org)./org.*100.0

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
# ╠═5ce58c17-ec49-4fa1-abed-961cb1380d6a
# ╠═78d8789f-952d-4abd-ba73-e540ad164e36
# ╠═31f05419-88f6-4dbb-9881-812864b1717f
# ╠═599cf329-ed37-453c-a190-b909847c1ec9
# ╠═7011dff6-2e07-4d03-b21f-94fd4359adc6
# ╠═9fbd5c15-1bb5-4344-b6a4-85549f1f8d7a
# ╠═1574adf8-9abe-49e3-8bfd-9970009a9373
# ╠═4f1b60f1-0a08-4983-a051-b278afbea4bb
# ╠═15f73a94-8d6f-4b1d-8314-076f8022f302
# ╠═b72343d0-e831-4699-9dcb-5251f4cf5738
# ╠═e20aeeed-5815-484c-acd7-ea57598b3642
# ╠═03f40097-f82b-4533-a3b4-9fbe3c3d67b3
# ╠═79d87deb-bf52-4d39-aae1-b16040d3ac62
# ╠═12cd49e8-9ed5-47e8-8415-5af39c49faf1
# ╠═829a7d93-e963-4084-b691-d9c12b41d0e2
# ╠═ad81e670-eec3-4610-92b7-99da8a9eb658
# ╠═3176cd9f-1c63-4b69-ab74-85a01a24f54a
# ╠═8cf3f13b-4c9d-4288-9be8-d53b0dfcb4ad
# ╠═d002c5fb-43b1-4873-8245-17fc94cde9c3
# ╠═42cbd777-f727-492d-b04f-49296f6be31f
# ╠═29ab9ff2-7472-4a95-80c8-6d64e4803a42
# ╠═ef392841-68b5-4407-8920-954cd4d78c4c
# ╠═71d799f2-b7e0-4145-9772-2fb13b0400a7
# ╠═672bd169-d9e7-4825-aa83-733d71e347d5
# ╠═6510d72a-01d7-46de-9385-e76b6167a8b2
# ╠═31229dd0-522d-4b49-8eb9-2a46255435cc
# ╠═199f10e6-2dfd-4fb1-9bf2-8b79f0023091
# ╠═58aa2d2d-cc91-4eea-bc45-3023c61f068a
# ╠═8eca087c-cb38-4c60-a574-f1161833af53
# ╠═49a8b09f-49e9-49fa-8aef-c9e7bb32f59c
# ╠═4161aef4-0b88-41dd-9f70-925989e2fa8a
# ╠═e4c126f3-e705-4f56-9443-026a708ccb90
# ╠═eedc9fdb-6c5a-441a-a437-e08864e2a202
# ╠═f948ecbc-f49f-40c7-b697-4a38cb95d458
# ╠═c7678500-8e16-46de-9583-4d36a0262304
# ╠═aff1eaa8-c680-407b-8a0f-dadae568db8b
# ╠═892cf99b-7672-4046-afdd-179673714235
# ╠═3784ab35-e63f-40a5-95ad-ff4183519ac6
# ╠═8aae0422-b0ec-4b61-88f2-1e7b7183a149
# ╠═7d0e83ca-4617-4d86-9266-507e0aa31d80
# ╠═8044d719-1fb1-4a50-87fc-45c5cf8de112
# ╠═884f9da5-0f32-423c-a71f-07abcd22f873
# ╠═6dab0dcf-c402-4bc8-a6f4-7c7cf305967f
# ╠═b380629a-880d-457d-bb2d-2cc6448ac4f5
# ╠═99d7ef10-1e01-4b0d-8bc3-98a761646149
# ╠═f68804fa-6f32-48cd-9df9-007302928a33
# ╠═07d0008d-06ee-42ba-8530-9440028d1746
# ╠═0af253f5-26ca-4423-84be-ab1f6165f1cf
# ╠═51f050ce-840e-4a06-af46-bfeb452dc30d
# ╠═abd86076-f11e-4171-9176-21c3779b3bc4
# ╠═2b509fc6-6ee9-44ca-bb95-df05246f2ca8
# ╠═c685bdc4-a9d0-4cbc-8f33-c7d560c10b5c
# ╠═10c316bc-6389-44dd-8a92-b772ffed9c84
# ╠═394deaba-a74e-45d0-9f02-c9ec7799e372
# ╠═4a909f87-f16b-4b88-bdc7-578e5849f01c
# ╠═26bc4084-4a58-4bbe-b0c2-4b087789a5c5
# ╠═dcc3a436-d079-47c1-808f-88179fa3d421
# ╠═7481f4bf-6622-495d-84de-d50d02b931f0
# ╠═fb68eaa0-ded3-40ac-90ab-ff0d21ed347e
# ╠═6d2408b0-2edf-48f5-aa13-83ff8d0c9955
# ╠═46751f67-d7ba-4e8b-8a82-c47becd9fcfd
# ╠═9f6122b4-ba52-4029-b86c-84320d7b6969
# ╠═f74ca574-f9d5-4b0f-b190-4a5a27ad78b1
# ╠═b7cae6f2-ffd7-4760-b3e6-1b869581d3b9
# ╠═f527f9c9-7bab-4fec-b0d0-2e24842a83af
# ╠═be703ad1-a05a-459f-a437-f02f97cef40b
# ╠═50a81a4a-a51f-4fb9-a104-07416fc479ae
# ╠═3d0f67c6-881e-4a77-b15f-be180f2bb567
# ╠═4e5b0aff-7c6f-46f9-8575-1d4a302973ef
# ╠═ab233b0f-62bd-45a0-aaef-d29bcdc7182e
# ╠═7d5bbbf8-6736-4e87-bae1-56e6fe4c4b22
# ╠═5a39f6a5-6c1a-4dbe-a6a7-15cbbca37ea8
# ╠═e4f421e0-516c-493e-ba62-36b93aaee21f
# ╠═6b255417-2e5e-4fe6-9acf-feda9c4713d0
# ╠═1c3ec334-4200-473e-a50e-8891f65f2110
# ╠═d7682fe5-6297-4220-965a-4a0a833e36b3
# ╠═9e923fe7-7059-4c4e-972c-36618ccc5149
# ╠═d9b5253c-96ac-4d45-a43e-8aabd3f54b49
# ╠═4fa6d902-8b70-4688-a3c2-f48e42919a66
# ╠═cd7a1778-0888-432a-b7ba-a5c0e0a94b32
# ╠═a2d60fc8-79c2-4c95-90b5-b83ebbf6357e
# ╠═e02f8b49-c372-40fe-80e6-af47f4a75d56
# ╠═c362bc14-4f0c-4111-b318-e640c96fb157
# ╠═6b0b29d4-276a-4888-8d17-a4468b98af95
