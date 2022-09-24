### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ bdaf89aa-42cd-4ed4-b009-a7820e233cc9
begin
	import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
	#include("/Users/janleppert/Documents/GitHub/ThermodynamicDataEstimator/src/ThermodynamicDataEstimator.jl")
	using ForwardDiff
	using GasChromatographySimulator
	using Plots
	using Optimization
	using OptimizationOptimJL
	using OptimizationBBO
	using ThermodynamicDataEstimator
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

!!! note

	This section and section "Simulate the test system" in separate notebook.
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
solutes = ["Octane", "2-Octanone", "2-Octanol", "1-Octanol"]

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

!!! note

	Make second notebook which uses real measurements as test data instead of simulations.
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
ii = 1

# ╔═╡ 322d274a-bfc7-46af-999a-6d3076f905c3
plot(GasChromatographySimulator.plot_temperature(par_meas[ii]), GasChromatographySimulator.plot_flow(par_meas[ii]))

# ╔═╡ 2446c0d3-8105-4644-a67d-ec39ec4dc1fc
begin
	pl_meas = Array{DataFrame}(undef, length(heating_rates))
	tR_meas = Array{Float64}(undef, length(heating_rates), length(solutes))
	tR_meas_randn = Array{Float64}(undef, length(heating_rates), length(solutes))
	TR_meas = Array{Float64}(undef, length(heating_rates), length(solutes))

	for i=1:length(heating_rates)
		pl_meas[i] = GasChromatographySimulator.simulate(par_meas[i])[1]
		for j=1:length(solutes)
			jj = findfirst(pl_meas[i].Name.==solutes[j])
			tR_meas[i,j] = pl_meas[i].tR[jj]
			tR_meas_randn[i,j] = pl_meas[i].tR[jj]*(1+randn()*0.005) # add here a random ± time 
			TR_meas[i,j] = pl_meas[i].TR[jj]
		end
	end
end

# ╔═╡ 5ce58c17-ec49-4fa1-abed-961cb1380d6a
tR_meas

# ╔═╡ 07a5cbf7-486d-4b6b-a346-6a11d01c21b0
tR_meas_randn

# ╔═╡ aee816b6-066c-43e0-a8aa-b301aea64543
TR_meas.+273.15

# ╔═╡ 9b69197d-d14d-4dd5-a1c5-f8e67b3e089e
md"""
## 1st estimates for parameters
Dimensionless heating rate for which ``T_{elu} ≈ T_{char}`` is ``r_T ≈ 0.6``. This heating rate is between the heating rates of the 4th (10°C/min) and 5th program (20°C/min). The mean value of the elution temperature at these two programs is used as a estimate for ``T_{char}``.

Based on the estimate for ``T_{char}`` an estimate for ``θ_{char}`` is calculated by [BlumbergBook]:

``
	θ_{char,e} = 22°C \left(\frac{T_{char,e}}{T_{st}}\right)^{0.7}
``

A good estimate for ``ΔCₚ`` is not yet available. For now it is set to a fixed value:

``
	ΔC_{p,e} = 100 J mol⁻¹ K⁻¹
``
"""

# ╔═╡ 84001173-491b-4cd9-9bde-1b3efb97a3c4
tMref = GasChromatographySimulator.holdup_time(150.0+273.15, 1.0/(60*1e6), 0.0, 30.0, 0.25e-3, "He", control="Flow")/60.0 # in min, approxymatly constant flow of 1mL/min

# ╔═╡ 75186010-825a-494c-b62a-1ce98284e6e5
dimless_heating_rates = heating_rates.*tMref./30.0

# ╔═╡ 430e4406-8d9c-42a5-99be-73893d2cbb65
heating_rates

# ╔═╡ 3cc1dc8c-04da-4f07-b2df-cf523c627dfa
begin
	Tchar_e = Array{Float64}(undef, size(TR_meas)[2])
	θchar_e = Array{Float64}(undef, size(TR_meas)[2])
	ΔCp_e = Array{Float64}(undef, size(TR_meas)[2])
	for i=1:size(TR_meas)[2]
		Tchar_e[i] = (TR_meas[4,i] + TR_meas[5,i])/2 + 273.15
		θchar_e[i] = 22.0*(Tchar_e[i]/273.15)^0.7
		ΔCp_e[i] = 100.0
	end
end

# ╔═╡ 8cef3271-02aa-4306-9352-499d835cd0f9
begin # true values of the Kcentric parameters
	Tchar = Array{Float64}(undef, length(solutes))
	θchar = Array{Float64}(undef, length(solutes))
	ΔCp = Array{Float64}(undef, length(solutes))
	for i=1:length(solutes)
		Tchar[i] = sub_def[i].Tchar
		θchar[i] = sub_def[i].θchar
		ΔCp[i] = sub_def[i].ΔCp
	end
end

# ╔═╡ da71a8b5-05c6-4911-8630-8b0363e4fbc7
begin
	plotly()
	plot_Tchar = scatter(1:size(TR_meas)[2], Tchar, label="true value", ylabel="Tchar in K")
	scatter!(plot_Tchar, 1:size(TR_meas)[2], Tchar_e, label="1st estimate")
	plot_θchar = scatter(1:size(TR_meas)[2], θchar, label="true value", ylabel="θchar in °C")
	scatter!(plot_θchar, 1:size(TR_meas)[2], θchar_e, label="1st estimate")
	plot_ΔCp = scatter(1:size(TR_meas)[2], ΔCp, label="true value", ylabel="ΔCp in °C")
	scatter!(plot_ΔCp, 1:size(TR_meas)[2], ΔCp_e, label="1st estimate")
	plot(plot_Tchar, plot_θchar, plot_ΔCp)
end

# ╔═╡ 599cf329-ed37-453c-a190-b909847c1ec9
md"""
## Define functions for optimization
"""

# ╔═╡ 2cf2e9ff-7222-4bfb-9fa4-1d565c9eb510
# optimize every solute separatly
function optimize_Kcentric_single(tR, L, d, gas, prog, opt, Tchar_e, θchar_e, ΔCp_e, methode, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp)
	optf = OptimizationFunction(ThermodynamicDataEstimator.opt_Kcentric, Optimization.AutoForwardDiff())
	opt_sol = Array{Any}(undef, size(tR)[2])
	for i=1:size(tR)[2]
		p = [tR[:,i], L, d, prog, opt, gas]
		x0 = [Tchar_e[i], θchar_e[i], ΔCp[i]]
		lb = [lb_Tchar[i], lb_θchar[i], lb_ΔCp[i]]
		ub = [ub_Tchar[i], ub_θchar[i], ub_ΔCp[i]]
		if methode == NelderMead() || methode == NewtonTrustRegion() || Symbol(methode) == Symbol(Newton())
			prob = OptimizationProblem(optf, x0, p)
		else
			prob = OptimizationProblem(optf, x0, p, lb=lb, ub=ub)
		end
		opt_sol[i] = solve(prob, methode)
	end
	return opt_sol
end

# ╔═╡ 286c8e0c-c4a2-4199-8a67-4da5f7da4380
methode=Optim.Newton()

# ╔═╡ 6d858656-7013-4423-8e0e-cc1441fdab02
Symbol(Optim.Newton()) === Symbol(Optim.Newton())

# ╔═╡ 8d05960e-3173-4b1f-9189-e08ee4634e90
Newton() == Optim.Newton()

# ╔═╡ 1b5ebfd3-a71b-437b-b65d-32ac59ab4610
(methode == NelderMead()) || (methode == NewtonTrustRegion()) || (Symbol(methode) == Symbol(Newton()))

# ╔═╡ fe9cca7f-d940-44b3-bdc4-0cecb0bc580e
methode == Optim.Newton()

# ╔═╡ 9120341e-3918-4bec-93ba-2ff048443b29
[(Tchar_e*0.9)[1], (θchar_e*0.8)[1], (ΔCp*0.5)[1]]

# ╔═╡ 8d84186c-b233-4c35-8982-7d6be6e36d43
opt = GasChromatographySimulator.Options(ng=true)

# ╔═╡ 572e036d-aefc-45ae-ab3f-6ce9a9cf3669
sol_single = optimize_Kcentric_single(tR_meas, 30.0, 0.25e-3, "He", prog_def, opt, Tchar_e, θchar_e, ΔCp_e, NelderMead(), Tchar_e*0.9, θchar_e*0.8, ΔCp_e*0.5, Tchar_e*1.1, θchar_e*1.2, ΔCp_e*1.5)

# ╔═╡ 09d7804f-93cb-4da7-a2f4-52a8561af6cf
# run different optimization methods in a script and save the results
# take some time to compute

# ╔═╡ 34f6518c-d720-4c24-a240-17ac63e81a84
#optimize_Kcentric_single(tR_meas, 30.0, 0.25e-3, "He", prog_def, opt, Tchar_e, θchar_e, ΔCp_e, LBFGS(), Tchar_e*0.9, θchar_e*0.8, ΔCp_e*0.5, Tchar_e*1.1, θchar_e*1.2, ΔCp_e*1.5)

# ╔═╡ 16a0a9ad-dc52-4458-87c6-fe82fb93ee01
#optimize_Kcentric_single(tR_meas, 30.0, 0.25e-3, "He", prog_def, opt, Tchar_e, θchar_e, ΔCp_e, BFGS(), Tchar_e*0.9, θchar_e*0.8, ΔCp_e*0.5, Tchar_e*1.1, θchar_e*1.2, ΔCp_e*1.5)

# ╔═╡ 7acc76d6-cc9f-41cb-9d3e-51e5570e284f
Tchar

# ╔═╡ 8dda11ba-574d-4e5f-a2bc-b6dc0b4c3fa4
θchar

# ╔═╡ 17d27066-b2d0-4e06-a86d-d9d36f89d907
ΔCp

# ╔═╡ 37061970-d234-4ad0-9ec2-39d13f0592e6
optimize_Kcentric_single(tR_meas, 30.0, 0.25e-3, "He", prog_def, opt, Tchar_e, θchar_e, ΔCp_e, NewtonTrustRegion(), Tchar_e*0.9, θchar_e*0.8, ΔCp_e*0.5, Tchar_e*1.1, θchar_e*1.2, ΔCp_e*1.5)

# ╔═╡ 492db22d-e4d2-4a89-9f4e-4fd4c7e68179
optimize_Kcentric_single(tR_meas, 30.0, 0.25e-3, "He", prog_def, opt, Tchar_e, θchar_e, ΔCp_e, Newton(), Tchar_e*0.9, θchar_e*0.8, ΔCp_e*0.5, Tchar_e*1.1, θchar_e*1.2, ΔCp_e*1.5)

# ╔═╡ 3c778f9c-1362-4d72-b920-5164333f0b9c
begin
	Tchar_opt = Array{Float64}(undef, length(solutes))
	θchar_opt = Array{Float64}(undef, length(solutes))
	ΔCp_opt = Array{Float64}(undef, length(solutes))
	for i=1:length(solutes)
		Tchar_opt[i] = sol_single[i][1]
		θchar_opt[i] = sol_single[i][2]
		ΔCp_opt[i] = sol_single[i][3]
	end
end

# ╔═╡ 495544fe-d6c9-4b5a-9662-05eb4ea2bebd
begin
	scatter!(plot_Tchar, 1:size(TR_meas)[2], Tchar_opt, label="optimized", m=:x, legend=:bottomright)
	scatter!(plot_θchar, 1:size(TR_meas)[2], θchar_opt, label="optimized", m=:x, legend=:bottomright)
	scatter!(plot_ΔCp, 1:size(TR_meas)[2], ΔCp_opt, label="optimized", m=:x, legend=:bottomright)
	plot(plot_Tchar, plot_θchar, plot_ΔCp)
end

# ╔═╡ 707bce9d-5cb1-40fe-8184-f2e6b2354910
# optimize all solutes together
function optimize_Kcentric_all(tR, L, d, gas, prog, opt, Tchar_e, θchar_e, ΔCp_e, methode, lb_Tchar, lb_θchar, lb_ΔCp, ub_Tchar, ub_θchar, ub_ΔCp)
	p = [tR, L, d, prog, opt, gas]
	x0 = [Tchar_e; θchar_e; ΔCp]
	lb = [lb_Tchar; lb_θchar; lb_ΔCp]
	ub = [ub_Tchar; ub_θchar; ub_ΔCp]
	optf = OptimizationFunction(ThermodynamicDataEstimator.opt_Kcentric, Optimization.AutoForwardDiff())
	if methode == NelderMead()
		prob = OptimizationProblem(optf, x0, p)
	else
		prob = OptimizationProblem(optf, x0, p, lb=lb, ub=ub)
	end
	opt_sol = solve(prob, methode)
	return opt_sol
end

# ╔═╡ 433ac08e-55d6-4298-b12c-f026201e2888
Tchar_e*(-Inf)

# ╔═╡ 88a01fcd-bbb7-4c81-8145-302f606912f0
Tchar_e*1.2

# ╔═╡ 62b70f05-6238-488c-b4c0-fa8043bcae32
sol_all = optimize_Kcentric_all(tR_meas, 30.0, 0.25e-3, "He", prog_def, opt, Tchar_e, θchar_e, ΔCp_e, NelderMead(), Tchar_e*(-Inf), θchar_e*(-Inf), ΔCp*(-Inf), Tchar_e*Inf, θchar_e*Inf, ΔCp*Inf)

# ╔═╡ 406b8abb-ae62-4324-8e78-86f0896adb0e
begin
	Tchar_opt_all = Array{Float64}(undef, length(solutes))
	θchar_opt_all = Array{Float64}(undef, length(solutes))
	ΔCp_opt_all = Array{Float64}(undef, length(solutes))
	for i=1:length(solutes)
		Tchar_opt_all[i] = sol_all[i]
		θchar_opt_all[i] = sol_all[4+i]
		ΔCp_opt_all[i] = sol_all[8+i]
	end
end

# ╔═╡ 9a68d8f3-444e-4d10-8098-a1cf18455d75
begin
	plotly()
	scatter!(plot_Tchar, 1:size(TR_meas)[2], Tchar_opt_all, label="optimized_all", m=:+, legend=:bottomright)
	scatter!(plot_θchar, 1:size(TR_meas)[2], θchar_opt_all, label="optimized_all", m=:+, legend=:bottomright)
	scatter!(plot_ΔCp, 1:size(TR_meas)[2], ΔCp_opt_all, label="optimized_all", m=:+, legend=:bottomright)
	plot(plot_Tchar, plot_θchar, plot_ΔCp)
end

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
# ╠═07a5cbf7-486d-4b6b-a346-6a11d01c21b0
# ╠═aee816b6-066c-43e0-a8aa-b301aea64543
# ╟─9b69197d-d14d-4dd5-a1c5-f8e67b3e089e
# ╠═84001173-491b-4cd9-9bde-1b3efb97a3c4
# ╠═75186010-825a-494c-b62a-1ce98284e6e5
# ╠═430e4406-8d9c-42a5-99be-73893d2cbb65
# ╠═3cc1dc8c-04da-4f07-b2df-cf523c627dfa
# ╠═8cef3271-02aa-4306-9352-499d835cd0f9
# ╠═da71a8b5-05c6-4911-8630-8b0363e4fbc7
# ╠═599cf329-ed37-453c-a190-b909847c1ec9
# ╠═2cf2e9ff-7222-4bfb-9fa4-1d565c9eb510
# ╠═286c8e0c-c4a2-4199-8a67-4da5f7da4380
# ╠═6d858656-7013-4423-8e0e-cc1441fdab02
# ╠═8d05960e-3173-4b1f-9189-e08ee4634e90
# ╠═1b5ebfd3-a71b-437b-b65d-32ac59ab4610
# ╠═fe9cca7f-d940-44b3-bdc4-0cecb0bc580e
# ╠═9120341e-3918-4bec-93ba-2ff048443b29
# ╠═8d84186c-b233-4c35-8982-7d6be6e36d43
# ╠═572e036d-aefc-45ae-ab3f-6ce9a9cf3669
# ╠═09d7804f-93cb-4da7-a2f4-52a8561af6cf
# ╠═34f6518c-d720-4c24-a240-17ac63e81a84
# ╠═16a0a9ad-dc52-4458-87c6-fe82fb93ee01
# ╠═7acc76d6-cc9f-41cb-9d3e-51e5570e284f
# ╠═8dda11ba-574d-4e5f-a2bc-b6dc0b4c3fa4
# ╠═17d27066-b2d0-4e06-a86d-d9d36f89d907
# ╠═37061970-d234-4ad0-9ec2-39d13f0592e6
# ╠═492db22d-e4d2-4a89-9f4e-4fd4c7e68179
# ╠═3c778f9c-1362-4d72-b920-5164333f0b9c
# ╠═495544fe-d6c9-4b5a-9662-05eb4ea2bebd
# ╠═707bce9d-5cb1-40fe-8184-f2e6b2354910
# ╠═433ac08e-55d6-4298-b12c-f026201e2888
# ╠═88a01fcd-bbb7-4c81-8145-302f606912f0
# ╠═62b70f05-6238-488c-b4c0-fa8043bcae32
# ╠═406b8abb-ae62-4324-8e78-86f0896adb0e
# ╠═9a68d8f3-444e-4d10-8098-a1cf18455d75
# ╠═6b0b29d4-276a-4888-8d17-a4468b98af95
