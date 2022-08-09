using GasChromatographySimulator

# read the input data
meas = DataFrame(CSV.File("data/measured_tR.csv", header=1))
n_solute = size(meas)[2]-10
# construct set of parameters
opt = GasChromatographySimulator.Options(OwrenZen5(), 1e-6, 1e-4, "inlet", true)
par_meas = Array{GasChromatographySimulator.Parameters}(undef, size(meas)[1])
for i=1:size(meas)[1]
    sys = GasChromatographySimulator.constructor_System(meas.L[i], meas.d[i], meas.df[i], meas.sp[i], meas.gas[i])
    prog = GasChromatographySimulator.constructor_Program(parse.(Float64, split(meas.timesteps[i])),
                                                            parse.(Float64, split(meas.tempsteps[i])),
                                                            parse.(Float64, split(meas.pinsteps[i])),
                                                            parse.(Float64, split(meas.poutsteps[i])),
                                                            sys.L)
    sub = Array{GasChromatographySimulator.Substance}(undef, n_solute)
    for j=1:n_solute
        # function for initial values of the thermodynamic parameters
        # Tchar_init -> calculate tM(@150°C) and based on initial φ₀ guess
        # determine which of the measured programs is close to the program,
        # where Telu ≈ Tchar (dimless heating rate 0.62) -> see
        # 'nb_heating_rate.jl'/'plnb_Telu_Tchar_sim_phase_ratio_corr.jl' in
        # VGGC-project, resp. in the ThermodynamicData-Project
        Tchar_init = 500.0
        θchar_init = 30.0
        ΔCp_init = 100.0
        sub[j] = GasChromatographySimulator.Substance(names(meas)[10+j],
                                                        "solute $(j)",
                                                        Tchar_init,
                                                        θchar_init,
                                                        ΔCp_init,
                                                        sys.a_df./sys.a_d,
                                                        "",
                                                        1e-4,
                                                        0.0,
                                                        0.0)
    end
    par_meas[i] = GasChromatographySimulator.Parameters(sys, prog, sub, opt)
end