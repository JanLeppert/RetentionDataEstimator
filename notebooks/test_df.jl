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
	using LeastSquaresOptim
	using PlutoUI
	TableOfContents()
end

# ╔═╡ 2c5e497e-09f2-4ffe-8df5-c455fd275771
md"""
# Test of Estimation of GC Parameters
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
			tR_meas_randn[i,j] = pl_meas[i].tR[jj]*(1+randn()*0.005) # add here a random ± time 
			tR_meas[i,j] = pl_meas[i].tR[jj]
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

# ╔═╡ b783d736-a4eb-447e-8f1f-636d5c497ee1
r_df(x, t, p, par_def, i) = GasChromatographySimulator.residency(x, t, par_def.prog.T_itp, par_def.prog.Fpin_itp, par_def.prog.pout_itp, par_def.col.L, par_def.col.d, p[1], par_def.col.gas, par_def.sub[i].Tchar, par_def.sub[i].θchar, par_def.sub[i].ΔCp, par_def.sub[i].φ₀; ng=par_def.opt.ng, vis=par_def.opt.vis, control=par_def.opt.control) # p[1] -> df

# ╔═╡ 367e0517-b188-4041-b4f4-6c7c6b7cbaba
function solving_migration_df(p, par_def, i)
    f_tx(t,p,x) = r_df(x, t, p, par_def, i)
    t₀ = 0.0
    xspan = (0.0, par_def.col.L)
    prob_tx = ODEProblem(f_tx, t₀, xspan, p)
    solution_tx = solve(prob_tx, alg=par_def.opt.alg, abstol=par_def.opt.abstol, reltol=par_def.opt.reltol)
	tR = solution_tx.u[end]
    return tR
end

# ╔═╡ 5f976e34-08c3-4292-adf0-73424aa0137b
function loss_df(tR, p, par_def)
    tR_p = Array{Float64}(undef, size(tR)[1], size(tR)[2])
    for i=1:size(tR)[1]
		for j=1:size(tR)[2]
        	tR_p[i,j] = solving_migration_df(p, par_def[i], j)
    	end
	end
    return sum((tR.-tR_p).^2), tR_p
end

# ╔═╡ 31229dd0-522d-4b49-8eb9-2a46255435cc
md"""
## Plot of the loss function
"""

# ╔═╡ d323849b-c2c9-416d-b195-2a1df161a1a9
retention_time = [measurements[:,2] measurements[:,3] measurements[:,5] measurements[:,4]]

# ╔═╡ b60a066e-7291-4722-bf9f-cc9e6499ed91
size(retention_time)

# ╔═╡ 13185458-7a35-43cc-b66f-5c7b76394147
retention_time_randn = [measurements_randn[:,2] measurements_randn[:,3] measurements_randn[:,5] measurements_randn[:,4]]

# ╔═╡ 1408a857-aa96-47b7-ba77-40d6ef7d5678
begin
	#L = 29.0:0.1:31.0
	#d = collect(0.24:0.002:0.26).*1e-3
	df = collect(0.1:0.01:0.3).*1e-6
	#retention_time = [measurements[:,2], measurements[:,3], measurements[:,5], measurements[:,4]] # attention to the correct order of retention times according to solutes? (or par.sub order)
	
	SSdf = Array{Float64}(undef, length(df))
	for k=1:length(df)
	    #for j=1:length(d)
	        #for i=1:length(L)
	            pp = [df[k]]
	            SSdf[k] = loss_df(retention_time, pp, par_meas)[1]
	        #end
	    #end
	end
end

# ╔═╡ aa47eca4-2ae4-406a-b308-bf1129e94718
begin
	plotly()
	plot(df, SSdf, xlabel="df in m")
end

# ╔═╡ 49a8b09f-49e9-49fa-8aef-c9e7bb32f59c
md"""
## Optimization Optim.jl

For now, just use the NelderMead() algorithm. 
"""

# ╔═╡ 78d8789f-952d-4abd-ba73-e540ad164e36
retention_time

# ╔═╡ df17fe3e-7ab8-46f1-929f-e17c037a6042
begin
	par_def = Array{GasChromatographySimulator.Parameters}(undef, length(heating_rates))
	for i=1:length(heating_rates)
		par_def[i] = GasChromatographySimulator.Parameters(col_def, prog_def[i], sub_def, opt_def)
	end
end

# ╔═╡ 4161aef4-0b88-41dd-9f70-925989e2fa8a
f_opt_df(x) = loss_df(retention_time, x, par_def)[1]

# ╔═╡ 8fa552fe-bd09-406d-a98e-188d4f4ce366
f_opt_df_randn(x) = loss_df(retention_time_randn, x, par_def)[1]

# ╔═╡ 155ff9ce-b6c0-408f-addb-25031832623e
p0_df = [0.23e-3]

# ╔═╡ 78a67009-b5bf-4dfd-8397-bf195bc00df0
res_opt_df_NM = Optim.optimize(f_opt_df, p0_df, NelderMead(), Optim.Options(iterations=5000))

# ╔═╡ a0521a79-cfc6-4c04-80e9-6bf61de778db
df_est = Optim.minimizer(res_opt_df_NM)

# ╔═╡ 8799a703-9d32-48a9-80dd-fe3fba7f3cae
res_opt_df_randn_NM = Optim.optimize(f_opt_df_randn, p0_df, NelderMead(), Optim.Options(iterations=5000))

# ╔═╡ 38827971-b433-4824-a3a6-6c0619201c95
df_est_randn = Optim.minimizer(res_opt_df_randn_NM)

# ╔═╡ 4e629dd7-d5da-4271-93ca-5d3da7105e14
#res_opt_df_LBFGS = Optim.optimize(f_opt_df, p0_df, LBFGS()) # runs forever

# ╔═╡ fad778df-9d71-427b-8842-3801cccc79ad
md"""
### Note

Optimization for film thickness with known Kcentric parameters (and known other GC parameters) is easily done with the NelderMead algorithm.
"""

# ╔═╡ 50e5f5f1-86b7-46ab-bdc7-18ad4b28bcf0
md"""
### Simulate with estimated ``d_f``
"""

# ╔═╡ 7ab0b7ee-1bca-4d44-89a0-f0f8a177bacc
begin
	par_est = Array{GasChromatographySimulator.Parameters}(undef, length(heating_rates))
	par_est_randn = Array{GasChromatographySimulator.Parameters}(undef, length(heating_rates))
	for i=1:length(heating_rates)
		col_est = GasChromatographySimulator.Column(col_def.L, col_def.d, df_est[1], col_def.sp, col_def.gas)
		col_est_randn = GasChromatographySimulator.Column(col_def.L, col_def.d, df_est_randn[1], col_def.sp, col_def.gas)
		par_est[i] = GasChromatographySimulator.Parameters(col_est, prog_def[i], sub_def, opt_def)
		par_est_randn[i] = GasChromatographySimulator.Parameters(col_est_randn, prog_def[i], sub_def, opt_def)
	end
end

# ╔═╡ 2e66c93c-3853-40f8-9648-b3785a70efe6
begin
	pl_est = Array{DataFrame}(undef, length(heating_rates))
	pl_est_randn = Array{DataFrame}(undef, length(heating_rates))
	tR_est = Array{Float64}(undef, length(heating_rates), length(solutes))
	tR_est_randn = Array{Float64}(undef, length(heating_rates), length(solutes))
	estimations = DataFrame(heating_rate=heating_rates)
	estimations_randn = DataFrame(heating_rate=heating_rates)
	for i=1:length(heating_rates)
		pl_est[i] = GasChromatographySimulator.simulate(par_est[i])[1]
		pl_est_randn[i] = GasChromatographySimulator.simulate(par_est_randn[i])[1]
		for j=1:length(solutes)
			jj = findfirst(pl_est[i].Name.==solutes[j])
			tR_est[i,j] = pl_est[i].tR[jj]
			tR_est_randn[i,j] = pl_est_randn[i].tR[jj]
		end
	end
	for j=1:length(solutes)
		estimations[!, Symbol(string("tR_", solutes[j]))] = tR_est[:,j]
		estimations_randn[!, Symbol(string("tR_", solutes[j]))] = tR_est_randn[:,j]
	end
	estimations, estimations_randn
end

# ╔═╡ 1b6952df-352e-415e-a08c-b256f7db9bda
measurements

# ╔═╡ a4608d23-ba51-4d24-aa4b-aa0c87f7a0c8
comparison = DataFrame(heating_rate=measurements.heating_rate, ΔtR_Octane=measurements[:,2]-estimations[:,2],ΔtR_2Octanone=measurements[:,3]-estimations[:,3],ΔtR_1Octanol=measurements[:,4]-estimations[:,4],ΔtR_2Octanol=measurements[:,5]-estimations[:,5])

# ╔═╡ 53a0208e-c464-4e07-b7aa-5fe5d1d8face
comparison_randn = DataFrame(heating_rate=measurements.heating_rate, ΔtR_Octane=measurements[:,2]-estimations_randn[:,2],ΔtR_2Octanone=measurements[:,3]-estimations_randn[:,3],ΔtR_1Octanol=measurements[:,4]-estimations_randn[:,4],ΔtR_2Octanol=measurements[:,5]-estimations_randn[:,5])

# ╔═╡ eedc9fdb-6c5a-441a-a437-e08864e2a202
md"""
## Gradient
"""

# ╔═╡ 0dda1a03-4b9b-4e71-bee1-4c9a34a31d4e
g(x) = ForwardDiff.gradient(f_opt_df, x)

# ╔═╡ cbeed1d3-310f-4717-9015-98ed6a9f092c
g([0.25])

# ╔═╡ 693159d1-aede-48ac-a660-7eb607f1ec3f
plot(0.1:0.01:0.4, )

# ╔═╡ 5733f75a-f2da-40b9-a0e9-4fa95e0ea5bb
#res_opt_Kcentric_BFGS = optimize(f_opt_Kcentric, p0_Kcentric, BFGS())

# ╔═╡ 779bf661-ca3e-465c-a7dc-24f611ed81bd
md"""
### Use IntervalRootFinding.jl on the gradient
"""

# ╔═╡ 5ece3d70-fcea-493d-9882-9695a16cdbe9
md"""
## Optimization LeastSquaresOptim.jl
"""

# ╔═╡ 2d7eb0fc-9084-4d1e-8aad-ee494af9a09a
p0_df

# ╔═╡ abc7fdd9-203a-4313-9aa3-dc6f65b6e0b2
res_LM = LeastSquaresOptim.optimize(f_opt_df, p0_df, LevenbergMarquardt())

# ╔═╡ 37ca7685-ad66-4a7f-9ed3-45a9cf1bbf57
res_Dog = LeastSquaresOptim.optimize(f_opt_df, p0_df, Dogleg())

# ╔═╡ df76a583-b639-40cf-ae07-e9547281c8ec
res_LM.minimizer # bad result

# ╔═╡ 7297a08c-ddb8-4bb2-8f63-7edddca64064
res_Dog.minimizer # better, but still worse than Optim with Nelder-Mead

# ╔═╡ f966b4ae-0e73-45b0-b37c-478a229c8b67
res_default = LeastSquaresOptim.optimize(f_opt_df, p0_df) # without definition of optimizer, it uses the default Optim.optimizer (here Nelder-Mead)

# ╔═╡ 431721c1-05cc-4346-af27-36310ac98286
res_default.minimizer

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
# ╠═b783d736-a4eb-447e-8f1f-636d5c497ee1
# ╠═367e0517-b188-4041-b4f4-6c7c6b7cbaba
# ╠═5f976e34-08c3-4292-adf0-73424aa0137b
# ╠═31229dd0-522d-4b49-8eb9-2a46255435cc
# ╠═b60a066e-7291-4722-bf9f-cc9e6499ed91
# ╠═d323849b-c2c9-416d-b195-2a1df161a1a9
# ╠═13185458-7a35-43cc-b66f-5c7b76394147
# ╠═1408a857-aa96-47b7-ba77-40d6ef7d5678
# ╠═aa47eca4-2ae4-406a-b308-bf1129e94718
# ╠═49a8b09f-49e9-49fa-8aef-c9e7bb32f59c
# ╠═78d8789f-952d-4abd-ba73-e540ad164e36
# ╠═df17fe3e-7ab8-46f1-929f-e17c037a6042
# ╠═4161aef4-0b88-41dd-9f70-925989e2fa8a
# ╠═8fa552fe-bd09-406d-a98e-188d4f4ce366
# ╠═155ff9ce-b6c0-408f-addb-25031832623e
# ╠═78a67009-b5bf-4dfd-8397-bf195bc00df0
# ╠═a0521a79-cfc6-4c04-80e9-6bf61de778db
# ╠═8799a703-9d32-48a9-80dd-fe3fba7f3cae
# ╠═38827971-b433-4824-a3a6-6c0619201c95
# ╠═4e629dd7-d5da-4271-93ca-5d3da7105e14
# ╠═fad778df-9d71-427b-8842-3801cccc79ad
# ╠═50e5f5f1-86b7-46ab-bdc7-18ad4b28bcf0
# ╠═7ab0b7ee-1bca-4d44-89a0-f0f8a177bacc
# ╠═2e66c93c-3853-40f8-9648-b3785a70efe6
# ╠═1b6952df-352e-415e-a08c-b256f7db9bda
# ╠═a4608d23-ba51-4d24-aa4b-aa0c87f7a0c8
# ╠═53a0208e-c464-4e07-b7aa-5fe5d1d8face
# ╠═eedc9fdb-6c5a-441a-a437-e08864e2a202
# ╠═0dda1a03-4b9b-4e71-bee1-4c9a34a31d4e
# ╠═cbeed1d3-310f-4717-9015-98ed6a9f092c
# ╠═693159d1-aede-48ac-a660-7eb607f1ec3f
# ╠═5733f75a-f2da-40b9-a0e9-4fa95e0ea5bb
# ╠═779bf661-ca3e-465c-a7dc-24f611ed81bd
# ╠═5ece3d70-fcea-493d-9882-9695a16cdbe9
# ╠═2d7eb0fc-9084-4d1e-8aad-ee494af9a09a
# ╠═abc7fdd9-203a-4313-9aa3-dc6f65b6e0b2
# ╠═37ca7685-ad66-4a7f-9ed3-45a9cf1bbf57
# ╠═df76a583-b639-40cf-ae07-e9547281c8ec
# ╠═7297a08c-ddb8-4bb2-8f63-7edddca64064
# ╠═f966b4ae-0e73-45b0-b37c-478a229c8b67
# ╠═431721c1-05cc-4346-af27-36310ac98286
# ╠═6b0b29d4-276a-4888-8d17-a4468b98af95
