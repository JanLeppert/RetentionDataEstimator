# Functions used to simulate chromatograms to use as test values for the optimization

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
    par_meas = Array{GasChromatographySimulator.Parameters}(undef, length(TPs))
    for i=1:length(TPs)
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