### A Pluto.jl notebook ###
# v0.19.10

using Markdown
using InteractiveUtils

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

# ╔═╡ 6afd1ed7-5e57-4ebd-8df3-894a6e208f41
prog_def # Array of the different programs

# ╔═╡ 2ae19586-4949-4610-9ff4-557ee54de1e8
opt_def

# ╔═╡ 1574adf8-9abe-49e3-8bfd-9970009a9373
opt_def_ = GasChromatographySimulator.Options(ng=true) # setting ng=true

# ╔═╡ 4f1b60f1-0a08-4983-a051-b278afbea4bb
function tR_calc(Tchar, θchar, ΔCp, L, d, df, prog, opt, gas)
	R = 8.3145
	k(x,t) = exp((ΔCp/R + Tchar/θchar)*(Tchar/prog.T_itp(x,t)-1) + ΔCp/R*log(prog.T_itp(x,t)/Tchar))
	rM(x,t) = GasChromatographySimulator.mobile_phase_residency(x,t, prog.T_itp, prog.Fpin_itp, prog.pout_itp, L, d, gas; ng=opt.ng, vis=opt.vis, control=opt.control)
	r(t,p,x) = (1+k(x,t))*rM(x,t)
	t₀ = 0.0
	xspan = (0.0, L)
	p = []
	prob = ODEProblem(r, t₀, xspan, p)
	solution = solve(prob, alg=opt.alg, abstol=opt.abstol, reltol=opt.reltol)
	tR = solution.u[end]
	return tR
end

# ╔═╡ e7cde70a-4530-4e50-a8c2-42b7ed9e86e1
col_def

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

# ╔═╡ 335e63bf-2679-4032-86b2-ef34e35927b7
df0 = 0.24e-6

# ╔═╡ 12cd49e8-9ed5-47e8-8415-5af39c49faf1
p0 = [Tchar0, θchar0, ΔCp0, L0, d0, df0]

# ╔═╡ 829a7d93-e963-4084-b691-d9c12b41d0e2
p0[1]

# ╔═╡ d002c5fb-43b1-4873-8245-17fc94cde9c3
function loss(tR, p, prog, opt)
	# loss function as sum over all programs and solutes
	gas = "He"
	Tchar0 = p[1] # Array length = number solutes
	θchar0 = p[2] # Array length = number solutes
	ΔCp0 = p[3] # Array length = number solutes
	L0 = p[4] # Float64
	d0 = p[5] # Float64
	df0 = p[6] # Float64
	tRcalc = Array{Float64}(undef, size(tR)[1], size(tR)[2])
	for j=1:size(tR)[2]
		for i=1:size(tR)[1]
			tRcalc[i,j] = tR_calc(Tchar0[j], θchar0[j], ΔCp0[j], L0, d0, df0, prog[i], opt, gas)
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
	df0 = p[6] # Float64
	tRcalc = Array{Float64}(undef, length(tR))
	for i=1:length(tR)
		tRcalc[i] = tR_calc(Tchar0, θchar0, ΔCp0, L0, d0, df0, prog[i], opt, gas)
	end
	return sum((tR.-tRcalc).^2), tRcalc
end	

# ╔═╡ 672bd169-d9e7-4825-aa83-733d71e347d5
p0_ = [Tchar0[1], θchar0[1], ΔCp0[1], L0, d0, df0]

# ╔═╡ 2290dd0b-a4fc-45d6-a5ad-55d43e827535
exp((ΔCp0[1]/8.3145 + Tchar0[1]/θchar0[1])*(Tchar0[1]/prog_def[1].T_itp(0.0,
	0.0)-1) + ΔCp0[1]*log(prog_def[1].T_itp(0.0,0.0)/Tchar0[1]))

# ╔═╡ 517ba539-fb70-43b7-8922-11e2092c11de
prog_def[1].T_itp(0.0,0.0)

# ╔═╡ bb107d6a-528a-47d2-ab1a-d6d0ef9aa422
tR_calc(Tchar0[1], θchar0[1], ΔCp0[1], L0, d0, df0, prog_def[7], opt_def_, "He")

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
## Compile optimization function
"""

# ╔═╡ 78d8789f-952d-4abd-ba73-e540ad164e36
retention_times = [measurements[:,2] measurements[:,3] measurements[:,4] measurements[:,5]]

# ╔═╡ 5a8219dc-182a-45f3-bdcb-276ebfaa035a
size(retention_times)

# ╔═╡ be41892d-710c-408c-9828-f2e06f2dcda0
tR_calc_ = Array{Float64}(undef, size(retention_times)[1], size(retention_times)[2])

# ╔═╡ 2755df8d-f552-4730-9472-4888102badef
tR_calc_

# ╔═╡ 9079658c-9a39-41f6-84b2-7ca8b028d551
for j=1:size(retention_times)[2]
		for i=1:size(retention_times)[1]
			tR_calc_[i,j] = tR_calc(Tchar0[j], θchar0[j], ΔCp0[j], L0, d0, df0, prog_def[i], opt_def_, "He")
		end
	end

# ╔═╡ ea8cb4f2-9cb7-4cc6-8f09-894c67108abb
collect(1:size(retention_times)[2])

# ╔═╡ 7159acb1-a756-4994-809e-270b5fce5b2b
sum((retention_times.-tR_calc_).^2)

# ╔═╡ 71d799f2-b7e0-4145-9772-2fb13b0400a7
loss(retention_times, p0, prog_def, opt_def)

# ╔═╡ b40cc773-2d65-4e73-b8ed-e99c5c2c4468
retention_times[:,1]

# ╔═╡ 6510d72a-01d7-46de-9385-e76b6167a8b2
loss_(retention_times[:,1], p0_, prog_def, opt_def)

# ╔═╡ 32d7f93f-414c-4541-8851-e30a33609afb
loss(retention_times, p0, prog_def, opt_def_)

# ╔═╡ 31f05419-88f6-4dbb-9881-812864b1717f
retention_times_randn = [measurements_randn[:,2] measurements_randn[:,3] measurements_randn[:,4] measurements_randn[:,5]]

# ╔═╡ 4161aef4-0b88-41dd-9f70-925989e2fa8a
f_opt_Kcentric(x) = loss_Kcentric(retention_times[:,1], x, col_def, prog_def, sub_def[1], opt_def)[1]

# ╔═╡ e4c126f3-e705-4f56-9443-026a708ccb90
f_opt_Kcentric_randn(x) = loss_Kcentric(retention_times_randn[:,1], x, col_def, prog_def, sub_def[1], opt_def)[1]

# ╔═╡ eedc9fdb-6c5a-441a-a437-e08864e2a202
md"""
## JuMP.jl
"""

# ╔═╡ 0c56fa69-5cdb-4c81-a663-b25d6c94e5cc
md"""
### Example

from: https://julia.quantecon.org/more\_julia/optimization\_solver\_packages.html
"""

# ╔═╡ f4a9bb1c-f281-4db1-9529-4d930657f5ab
# solve
# max( x[1] + x[2] )
# st sqrt(x[1]^2 + x[2]^2) <= 1

# ╔═╡ 465ef586-06cf-4fde-88b0-2320527f022a
function squareroot(x) # pretending we don't know sqrt()
    z = x # Initial starting point for Newton’s method
    while abs(z*z - x) > 1e-13
        z = z - (z*z-x)/(2z)
    end
    return z
end

# ╔═╡ fdc5afae-0d6b-4ced-b4c8-6a394a2e56a7
m = Model(Ipopt.Optimizer)

# ╔═╡ 6c4babe7-fc1c-4578-8022-7b24b673e779
# need to register user defined functions for AD
JuMP.register(m,:squareroot, 1, squareroot, autodiff=true)

# ╔═╡ 799018d1-f85a-473f-a71b-396e995a7d79
@variable(m, x[1:2], start=0.5) # start is the initial condition

# ╔═╡ 690f3dce-e74e-4e4e-bd17-e2a3bcfe37c6
@objective(m, Max, sum(x))

# ╔═╡ 8e74230e-c0cf-4958-9e5a-182a3cd59f28
@NLconstraint(m, squareroot(x[1]^2+x[2]^2) <= 1)

# ╔═╡ b56eaa46-fa4d-43ef-89de-005f5490021e
@show JuMP.optimize!(m)

# ╔═╡ c7678500-8e16-46de-9583-4d36a0262304
md"""
### My optimization
"""

# ╔═╡ 12bb4759-31ef-4ba8-8b0e-4b83dd56a1ed
model = Model(Ipopt.Optimizer)

# ╔═╡ bec272ed-ab3b-4f75-9b1b-4b42626808f7
@variable(model, 300 <= xchar <= 400)

# ╔═╡ b40c9046-6cf3-4898-a6fd-cdea92eed641
@variable(model, 20 <= ychar <= 50)

# ╔═╡ c7d74ae3-0c00-4105-a6bd-61ea84e79cce
@variable(model, 0 <= Δxp <= 200)

# ╔═╡ 646903d4-62e7-4657-8b31-499ecc4cdbbf
f_opt(x,y,z) = loss_(retention_times[:,1], [x,y,z,30.0,0.25e-3,0.25e-6], prog_def, opt_def_)[1]

# ╔═╡ 3edfab78-c28c-4afe-b1f2-2c0443031585
f_opt(Tchar0[1], θchar0[1], ΔCp0[1]) 

# ╔═╡ f5aff786-f2ea-423e-81d5-35ecac692b96
f_opt_(xyz) = f_opt(xyz[1],xyz[2],xyz[3])

# ╔═╡ 76717785-0f1d-4c2c-9250-053ac4d309fb
grad(xyz) = ForwardDiff.gradient(f_opt_, xyz)

# ╔═╡ 78cb372a-85f4-4e56-bde4-6ae55c70c91d
rM(xt) = GasChromatographySimulator.mobile_phase_residency(xt[1], xt[2], prog_def[1].T_itp, prog_def[1].Fpin_itp, prog_def[1].pout_itp, L0, d0, "He"; ng=opt_def_.ng, vis=opt_def_.vis, control=opt_def_.control)

# ╔═╡ 7c8cd3ee-8898-4c81-98a3-bd9c00672b9b
∂rM(xt) = ForwardDiff.gradient(rM, xt)

# ╔═╡ fedf0444-feab-41fa-8db6-b851cc190f5d
∂rM([0.0, 0.0])

# ╔═╡ c48f6289-9c19-4fa4-b525-10af808af240
grad([Tchar0[1], θchar0[1], ΔCp0[1]])

# ╔═╡ bf5d74ee-0d83-43a1-b38b-4ce452139eb5
g_opt_Kcentric(x,y,z) = f_opt_Kcentric([x, y, z])

# ╔═╡ 02dcc8bd-9bee-4f84-a3e9-b20648f2ec2b
register(model, :f_opt, 3, f_opt; autodiff=true)

# ╔═╡ f9f8bb48-27ba-467c-b7f2-c7d8fc5ecc26
@NLobjective(model, Min, f_opt(xchar, ychar, Δxp))

# ╔═╡ 3d2579de-ef06-4f68-8afb-43ffae030260
print(model)

# ╔═╡ ff4980b1-fb1d-4e4d-b9ef-3b1d12d0a586
#optimize!(model)

# ╔═╡ 6b0b29d4-276a-4888-8d17-a4468b98af95
md"""
# End
"""

# ╔═╡ bdaf89aa-42cd-4ed4-b009-a7820e233cc9
begin
	import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
	include("/Users/janleppert/Documents/GitHub/ThermodynamicDataEstimator/src/ThermodynamicDataEstimator.jl")
	using ForwardDiff
	using GasChromatographySimulator
	using Plots
	using JuMP
    using Ipopt
	using PlutoUI
	TableOfContents()
end

# ╔═╡ daf8edd9-b2b5-4def-8102-a3428410ce9f
using ForwardDiff

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
# ╠═6afd1ed7-5e57-4ebd-8df3-894a6e208f41
# ╠═2ae19586-4949-4610-9ff4-557ee54de1e8
# ╠═1574adf8-9abe-49e3-8bfd-9970009a9373
# ╠═4f1b60f1-0a08-4983-a051-b278afbea4bb
# ╠═e7cde70a-4530-4e50-a8c2-42b7ed9e86e1
# ╠═15f73a94-8d6f-4b1d-8314-076f8022f302
# ╠═b72343d0-e831-4699-9dcb-5251f4cf5738
# ╠═e20aeeed-5815-484c-acd7-ea57598b3642
# ╠═03f40097-f82b-4533-a3b4-9fbe3c3d67b3
# ╠═79d87deb-bf52-4d39-aae1-b16040d3ac62
# ╠═335e63bf-2679-4032-86b2-ef34e35927b7
# ╠═12cd49e8-9ed5-47e8-8415-5af39c49faf1
# ╠═829a7d93-e963-4084-b691-d9c12b41d0e2
# ╠═5a8219dc-182a-45f3-bdcb-276ebfaa035a
# ╠═d002c5fb-43b1-4873-8245-17fc94cde9c3
# ╠═ef392841-68b5-4407-8920-954cd4d78c4c
# ╠═be41892d-710c-408c-9828-f2e06f2dcda0
# ╠═9079658c-9a39-41f6-84b2-7ca8b028d551
# ╠═2755df8d-f552-4730-9472-4888102badef
# ╠═ea8cb4f2-9cb7-4cc6-8f09-894c67108abb
# ╠═7159acb1-a756-4994-809e-270b5fce5b2b
# ╠═71d799f2-b7e0-4145-9772-2fb13b0400a7
# ╠═b40cc773-2d65-4e73-b8ed-e99c5c2c4468
# ╠═672bd169-d9e7-4825-aa83-733d71e347d5
# ╠═6510d72a-01d7-46de-9385-e76b6167a8b2
# ╠═32d7f93f-414c-4541-8851-e30a33609afb
# ╠═2290dd0b-a4fc-45d6-a5ad-55d43e827535
# ╠═517ba539-fb70-43b7-8922-11e2092c11de
# ╠═bb107d6a-528a-47d2-ab1a-d6d0ef9aa422
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
# ╠═eedc9fdb-6c5a-441a-a437-e08864e2a202
# ╠═0c56fa69-5cdb-4c81-a663-b25d6c94e5cc
# ╠═f4a9bb1c-f281-4db1-9529-4d930657f5ab
# ╠═465ef586-06cf-4fde-88b0-2320527f022a
# ╠═fdc5afae-0d6b-4ced-b4c8-6a394a2e56a7
# ╠═6c4babe7-fc1c-4578-8022-7b24b673e779
# ╠═799018d1-f85a-473f-a71b-396e995a7d79
# ╠═690f3dce-e74e-4e4e-bd17-e2a3bcfe37c6
# ╠═8e74230e-c0cf-4958-9e5a-182a3cd59f28
# ╠═b56eaa46-fa4d-43ef-89de-005f5490021e
# ╠═c7678500-8e16-46de-9583-4d36a0262304
# ╠═12bb4759-31ef-4ba8-8b0e-4b83dd56a1ed
# ╠═bec272ed-ab3b-4f75-9b1b-4b42626808f7
# ╠═b40c9046-6cf3-4898-a6fd-cdea92eed641
# ╠═c7d74ae3-0c00-4105-a6bd-61ea84e79cce
# ╠═646903d4-62e7-4657-8b31-499ecc4cdbbf
# ╠═3edfab78-c28c-4afe-b1f2-2c0443031585
# ╠═daf8edd9-b2b5-4def-8102-a3428410ce9f
# ╠═f5aff786-f2ea-423e-81d5-35ecac692b96
# ╠═76717785-0f1d-4c2c-9250-053ac4d309fb
# ╠═78cb372a-85f4-4e56-bde4-6ae55c70c91d
# ╠═7c8cd3ee-8898-4c81-98a3-bd9c00672b9b
# ╠═fedf0444-feab-41fa-8db6-b851cc190f5d
# ╠═c48f6289-9c19-4fa4-b525-10af808af240
# ╠═bf5d74ee-0d83-43a1-b38b-4ce452139eb5
# ╠═02dcc8bd-9bee-4f84-a3e9-b20648f2ec2b
# ╠═f9f8bb48-27ba-467c-b7f2-c7d8fc5ecc26
# ╠═3d2579de-ef06-4f68-8afb-43ffae030260
# ╠═ff4980b1-fb1d-4e4d-b9ef-3b1d12d0a586
# ╠═6b0b29d4-276a-4888-8d17-a4468b98af95
