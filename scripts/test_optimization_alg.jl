using ForwardDiff
using GasChromatographySimulator
using Plots
using Optimization
using OptimizationOptimJL
using OptimizationBBO
using Dierckx

# test system (put this in a function)
# function to create a conventional GC-system to calculate test chromatograms
"""
    conventional_GC(L, d, df, sp, gas, TP, PP, solutes, db_path, db_file)

Description.
"""
function conventional_GC(L, d, df, sp, gas, TP, PP, solutes, db_path, db_file)
	opt = GasChromatographySimulator.Options()
	col = GasChromatographySimulator.Column(L, d, df, sp, gas)
	prog = GasChromatographySimulator.Program(TP, PP, L; pout="vacuum", time_unit="min")
	sub = GasChromatographySimulator.load_solute_database(db_path, db_file, sp, gas, solutes, zeros(length(solutes)), zeros(length(solutes)))
	par = GasChromatographySimulator.Parameters(col, prog, sub, opt)
	return par
end

"""
    sim_test_chrom(L, d, df, sp, gas, TPs, PPs, solutes, db_path, db_file)

Description.
"""
function sim_test_chrom(L, d, df, sp, gas, TPs, PPs, solutes, db_path, db_file)
    par_meas = Array{GasChromatographySimulator.Parameters}(undef, length(heating_rates))
    for i=1:length(heating_rates)
        par_meas[i] = conventional_GC(L, d, df, sp, gas, TPs[i], PPs[i], solutes, db_path, db_file)
    end

    #pl_meas = Array{DataFrame}(undef, length(TPs))
    tR_meas = Array{Float64}(undef, length(TPs), length(solutes))
    tR_meas_randn = Array{Float64}(undef, length(TPs), length(solutes))
    for i=1:length(TPs)
        pl_meas = GasChromatographySimulator.simulate(par_meas[i])[1]
        for j=1:length(solutes)
            jj = findfirst(pl_meas.Name.==solutes[j])
            tR_meas[i,j] = pl_meas.tR[jj]
            tR_meas_randn[i,j] = pl_meas.tR[jj]*(1+randn()*0.005) # add here a random ± time
        end
    end
    return tR_meas, tR_meas_randn, par_meas
end

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
tR_meas, tR_meas_randn, par_meas = sim_test_chrom(L, d, df, sp, gas, TPs, PPs, solutes, db_path, db_file)


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

"""
    reference_holdup_time(prog, L, d, fas; control="Pressure")

Description.
"""
function reference_holdup_time(prog, L, d, gas; control="Pressure")
    Tref = 150.0
    # estimate the time of the temperature program for T=Tref
    t_ = prog.time_steps
    T_ = prog.temp_steps
    spl = Spline1D(cumsum(t_), T_ .- Tref)
    tref = roots(spl)[1]
    # inlet and outlet pressure at time tref
    Fpin_ref = prog.Fpin_itp(tref)
    pout_ref = prog.pout_itp(tref)
    # hold-up time calculated for the time of the program, when T=Tref
    tMref = GasChromatographySimulator.holdup_time(Tref+273.15, Fpin_ref, pout_ref, L, d, gas, control=control)
    return tMref
end

tMref_def = Array{Float64}(undef, length(TPs))
for i=1:length(TPs)
    tMref_def[i] = reference_holdup_time(par_def[i].prog, L_def, d_def, gas; control="Pressure")/60.0
end 
# 1. Start values of the parameters