using CSV, DataFrames, Statistics

#include("Inomura2020.jl");

function generate_params(PIcurve)
    ESDCell = 0.7   # Prochlorococcus ESD ~0.5-0.7 um
    RadCell = ESDCell / 2.0 
    # volume of cell (um^3)
    VolCell = (4.0*3.14159/3.0)*RadCell^3.0
    # surface area of cell (um^2)
    SurfAreaCell = (4.0*3.14159)*RadCell^2.0
    # Alternatively, carbon quota for Prochlorococcus from direct measurement
    Q_C = 0.0609/12.0  # pmol C cell-1,  MED-4, P-limited, Bertilsson et al (2003)

    # Use Redfield ratio to evaluate Q_N and Q_P
    Q_N = (16.0/106.0)*Q_C           # pmol N cell-1
    Q_P = (1.0/106.0)*Q_C            # pmol P cell-1
    # Hawco et al report 1.6*10^-5 to 4.0*10^5 Fe:C ratios in Prochlorococcus
    Q_Fe = 1.6e-5*Q_C                # pmol Fe cell-1

    # set allometric scaling coeffs for Vmax_NO3 and K_NO3 (Edwards et al, 2012)
    a_V_NO3  = 10.0^(-8.1)      # umol N cell-1 day-1
    b_V_NO3 = 0.82              # power
    a_K_NO3 = 10.0^(-0.84)      # umol N L-1
    b_K_NO3 = 0.33              # power
    Vmax_NO3 = a_V_NO3 * (VolCell^b_V_NO3)    # umol N cell-1 day-1
    K_NO3 = a_K_NO3 * (VolCell^b_K_NO3)       # umol N L-1
    # evaluate N specific max N uptake
    Vmax_NO3_specific = Vmax_NO3 / (Q_N * 1.0e-6)    # day-1

    # set allometric scaling coeffs for Vmax_PO4 and K_PO4 (Edwards et al, 2012)
    a_V_PO4  = 10.0^(-9.1)      # umol N cell-1 day-1
    b_V_PO4 = 1.0               # power
    a_K_PO4 = 10.0^(-1.4)      # umol N L-1
    b_K_PO4 = 0.41              # power
    Vmax_PO4 = a_V_NO3 * (VolCell^b_V_PO4)    # umol N cell-1 day-1
    K_PO4 = a_K_PO4 * (VolCell^b_K_PO4)       # umol N L-1
    # evaluate N specific max N uptake
    Vmax_PO4_specific = Vmax_PO4 / (Q_P * 1.0e-6)    # day-1

    # set allometric scaling for Fe uptake, Shaked et al iron uptake, (Lis et al, Life, 2015; Lis et al, ISMEj 2015)
    Kin_SA_Fe = 10.0e-10                  # L day-1 um-2 
    Kin_Fe = Kin_SA_Fe * SurfAreaCell    # L cell-1 day-1
    Q_Fe_molar = Q_Fe * 1.0e-12 # Q_Fe_molar - convert to pmol to mol Fe cell-1
    Kin_Fe_Specific = Kin_Fe / (Q_Fe_molar)  # L cell-1 day-1 (mol Fe cell-1)-1 = day-1 / (mol Fe L-1) 
    
    mort = 0.8 # constant death rate
    K_R = 0.05 # respiration rate, day-1
    
    # Partensky et al 1993 .............................
    alpha_Chl_Partensky  = 0.041      # trying out some representative values
    # I_b and I_m are acclimated characteristics... depend on recent light
    # reported by Partensky
    I_m = 341.0             # 
    I_b_Partensky = 730.0             # I_b = P_S_Chl / beta_Chl 
    # Partensky don't report P_S_Chl... estimate to roughly match data
    P_S_Chl_Partensky = 15.0    # fg C (fg Chl)-1 hour-1
    beta_Chl_Partensky  = P_S_Chl_Partensky / I_b_Partensky
    # ....................................................

    # Moore and Chisholm (1999) ...........................
    # 9211 Low light strain, adapted at 9 µmol quanta m-2 s-1
    alpha_Chl_MC_9211_Lo = 0.09   # estimated directly from bar graph
    P_S_Chl_MC_9211_Lo = 3.0      # guesstimated
    beta_Chl_MC_9211_Lo = 0.010   # guesstimated
    P_S_Chl_Lo_MM = 97.47         # mol C (mol C in Chl)-1 day-1
    # 9215 High light strain, adapted at 70 µmol quanta m-2 s-1
    alpha_Chl_MC_9215_Hi = 0.05   # estimated directly from bar graph
    P_S_Chl_MC_9215_Hi = 8.0      # guesstimated
    beta_Chl_MC_9215_Hi = 0.008   # guesstimated
    P_S_Chl_Hi_MM = 259.92        # mol C (mol C in Chl)-1 day-1
    
     if PIcurve == "Pa"
        alpha_Chl = alpha_Chl_Partensky
        P_S_Chl = P_S_Chl_Partensky
        beta_Chl = beta_Chl_Partensky
    end
    if PIcurve == "LL"
        alpha_Chl = alpha_Chl_MC_9211_Lo 
        P_S_Chl = P_S_Chl_MC_9211_Lo
        beta_Chl = beta_Chl_MC_9211_Lo
        nu_I_max = P_S_Chl_Lo_MM
    end
    if PIcurve == "HL"
        alpha_Chl = alpha_Chl_MC_9215_Hi 
        P_S_Chl = P_S_Chl_MC_9215_Hi
        beta_Chl = beta_Chl_MC_9215_Hi
        nu_I_max = P_S_Chl_Hi_MM
    end
    # set I_b accordingly 
    I_b = P_S_Chl / beta_Chl 
    
    # MacroMolecule Allocation parameters from Inomura et al. 2020
    A_Pho = 1.60e1               # mol C (mol C in Chl)-1
    A_Bio = 2.71e-1              # mol C (mol C)-1 day-1
    A_RNA_P = 4.23e-3            # mol P (mol C)-1 day
    A_Pho_PChl = 2.83e-2         # mol P (mol C in Chl)-1
    Q_C_Pro_Other = 2.40e-1     # mol C (mol C)-1
    Q_C_Other0 = 1.82e-2         # mol C (mol C)-1
    Q_C_DNA = 9.41e-4            # mol C (mol C)-1
    Q_Pmin_RNA = 2.23e-4         # mol P (mol C)-1
    Y_RNA_CP = 10.7/1.0
    Y_Plip_CP = 40.0/1.0
    
    A_I = alpha_Chl/P_S_Chl
    E = 7.74e-1                  # synth cost dimensionless
    m = 3.93e-1 #K_R             # maintenance respiration rate, 3.93e-1 in Inomura et al.
    
    return (VNmax = Vmax_NO3_specific, VPmax= Vmax_PO4_specific, 
        KinFe = Kin_Fe_Specific, KN = K_NO3, KP = K_PO4, KR = K_R, mort = mort,
        alphaChl = alpha_Chl, PSchl = P_S_Chl, Ib = I_b, Qc = Q_C, 
        A_Pho = A_Pho, A_Bio = A_Bio, A_RNA_P = A_RNA_P, A_Pho_PChl = A_Pho_PChl,
        Q_C_Pro_Other = Q_C_Pro_Other, Q_C_Other0 = Q_C_Other0, Q_C_DNA = Q_C_DNA,
        Q_Pmin_RNA = Q_Pmin_RNA, Y_RNA_CP = Y_RNA_CP, Y_Plip_CP = Y_Plip_CP, 
        A_I = A_I, nu_I_max = nu_I_max, E = E, m = m)
end

function nut_min_uptake(input, params)
    VN = params.VNmax .* input.DIN ./ (input.DIN .+ params.KN)
    VP = params.VPmax .* input.PO4 ./ (input.PO4 .+ params.KP)
    VFe= params.KinFe .* input.Fe_d .* 1.0e-9  # because [Fe_d] in nmoles L-1
    
    return min(VN, VFe) # per day
end

function P_Chl(input, params, ::Nothing; photo_inhib = true)
    P_Chl = params.PSchl .* (1.0 .- exp.(-params.alphaChl .* input.PARav ./ params.PSchl))
    if photo_inhib
        P_Chl .*= exp.(-1.0 .* input.PARav ./ params.Ib)
    end
    return P_Chl
end

function P_Chl(input, params, hPARfrac::Array; photo_inhib = true)
    nhour = length(hPARfrac)
    nz = size(input,1)
    I_hourly = zeros(nz, nhour)
    for k in 1:nhour
        I_hourly[:,k] = input.PARmax.*  hPARfrac[k]
    end
    
    P_Chl_hourly = params.PSchl .* (1.0 .- exp.(-params.alphaChl .* I_hourly ./ params.PSchl))
    if photo_inhib
        P_Chl_hourly .*= exp.(-1.0 .* I_hourly ./ params.Ib)
    end
    
    P_Chl = mean(P_Chl_hourly, dims=2)[:,1]
    return P_Chl # per hour
end

# nu_I(I), photosynthesis vs light. Evaluate A_Chl and B_Chl
function photosynth(I,pa)
    nu_I = pa.nu_I_max*(1.0 - exp(-pa.A_I*I))
    A_Chl = (1.0 + pa.E)/nu_I
    B_Chl = pa.m / nu_I
    return nu_I, A_Chl, B_Chl
end 

# mu_max_I, max light-limited growth rate for given I (day-1)
function mu_max_I_eval(A_Chl,B_Chl,pa)
    # coefficients of quadratic
    a_M = pa.Y_RNA_CP*pa.A_RNA_P*(pa.A_Pho*A_Chl + pa.A_Bio)
    b_M = (  (1.0 + pa.A_Pho + pa.Y_Plip_CP*pa.A_Pho_PChl)*A_Chl + pa.A_Bio
                  + pa.Y_RNA_CP*pa.A_RNA_P*(pa.A_Pho*B_Chl + pa.Q_C_Pro_Other)  )
    Q_C_Other = pa.Q_C_Other0 + pa.Q_C_DNA + pa.Q_C_Pro_Other
    c_M = ( (1.0 + pa.A_Pho + pa.Y_Plip_CP*pa.A_Pho_PChl)*B_Chl
                  + Q_C_Other + pa.Y_RNA_CP*pa.Q_Pmin_RNA - 1.0     )
    # evaluate root of quadratic
    mu_max_I = (- b_M + (b_M*b_M - 4.0*a_M*c_M)^0.5 )/(2.0*a_M)
    return mu_max_I
end

# Chlorophyll to carbon ratio (moles C in Chl (mol C)-1)
function ChlC(mu,A_Chl,B_Chl)
    Q_C_Chl =  A_Chl*mu + B_Chl
    return Q_C_Chl
end

function calc_mu(input, params; hPARfrac = nothing, photo_inhib = true)
    nz = size(input, 1)
    
    A_Chl_out  = zeros(nz)
    B_Chl_out  = zeros(nz)
    mumaxI_out = zeros(nz)
    CN         = zeros(nz)
    CP         = zeros(nz)
    Chl_C      = zeros(nz)
    mu_C       = zeros(nz)
    mu_nut     = zeros(nz)
    
    # Factor to convert Chl-a:C from mol C in Chl-a (mol C)-1 [Inomura]  to  fg Chl-a (fg C)-1 [Partensky]
    Mchl=893.49  # molar mass of chlorophyll (g mol-1)
    Mc=12.0      # molar mass of carbon (g mol-1)
    C_per_Chla = 55.0   # Chlorophyll-a =  C₅₅H₇₂O₅N₄Mg
    Chl_factor = Mchl /  (Mc * C_per_Chla)  
    
    # growth-rate INDEPENDENT photosynth params: =================
    # Use DAILY MEAN light to evaluate these... 
    # ASSUMES ACCLIMATION HAPPENS ON DAILY TIMESCALE
    # otherwise acclimates Chl:C on hourly scale which leads to odd
    for i in 1:nz
        A_Chl_out[i] = photosynth(input.PARav[i],params)[2]
        B_Chl_out[i] = photosynth(input.PARav[i],params)[3]
    end
    
    # Iterative determination of Chl:C, growth rate, biomass
    mu_auto = (1.5 .*  exp.(-input.Depth/150.0))   # first "guess" for mu for initial Chl:C
    for j in 1:40
        for i in 1:nz
            mumaxI_out[i] = mu_max_I_eval(A_Chl_out[i], B_Chl_out[i], params)
            Chl_C[i] = ChlC(mu_auto[i],A_Chl_out[i],B_Chl_out[i]) * Chl_factor
            #CN[i] = 1.0 / NC_eval(mu_auto[i],A_Chl_out[i],B_Chl_out[i],pa)
            #CP[i] = 1.0 / PC_eval(mu_auto[i],A_Chl_out[i],B_Chl_out[i],pa)
        end
        
        # Eval photosynth (d-1): fg C (fg Chl)-1 h-1 * fg Chl (fg C)-1 * h d-1
        Pcalc = P_Chl(input, params, hPARfrac; photo_inhib = photo_inhib) .* Chl_C .* 24.0  # estimated photosynth growth 
        mumax = max.(mumaxI_out, 0.0)  # Max possible light-lim growth rate
        #P = min.(Pcalc, Pmax)         # Photosynth "growth rate" 
        
        # maintenance respiration rate used here and in Inomura et al are quite different, almost 10 times
        mu_C .= min.(mumax, Pcalc .- params.KR)         # net carbon growth rate
        mu_nut .= nut_min_uptake(input, params) # nutrient limited growth rate
        mu_min = min.(mu_C , mu_nut)  # minimal growth rate
        
        mu_auto .= max.(mu_min, 0,0)
    end
    
    # evaluate heterotrophic growth rate by the difference of observed growth rate and autotrophic growth rate
    mu_het = input.Pro_mu .- mu_auto
    mu_het = max.(mu_het, 0.0)
    
    # evaluate biomass from observed growth rate and quadratic grazing
    B_tot = input.Pro_mu ./ params.mort 
    B_Chl = B_tot .* Chl_C * 12.0   # Chl concentration,  ug Chl L-1
    # estimate cell density # /10^5 cells ml-1
    CD_tot = 1.0e-5 * 1.0e3 .* B_tot ./ params.Qc   # 1000.0* micromol C L-1 / picomol C cell-1
    
    # evaluate biomass from simulated autotrophic growth rate and quadratic grazing
    B_auto = mu_auto ./ params.mort 
    CD_auto = 1.0e-5 * 1.0e3 .* B_auto ./ params.Qc   # 1000.0* micromol C L-1 / picomol C cell-1
    
    return (mu_auto = mu_auto, B = B_tot, Chl = B_Chl, mu_C = mu_C, mu_nut = mu_nut, CD_tot = CD_tot, mu_het = mu_het, CD_auto = CD_auto)
end

function output_summary(input, hPARfrac)
    params_HL = generate_params("HL")
    params_LL = generate_params("LL")
    
    out1 = calc_mu(input, params_HL, hPARfrac = hPARfrac, photo_inhib = true)
    out2 = calc_mu(input, params_HL, hPARfrac = hPARfrac, photo_inhib = false)
    out3 = calc_mu(input, params_LL, hPARfrac = hPARfrac, photo_inhib = true)
    out4 = calc_mu(input, params_LL, hPARfrac = hPARfrac, photo_inhib = false)
    
    out = (out1, out2, out3, out4)
    
    nz = length(out1.mu_auto)
    
    output = zeros(nz, 4, 5)

    for i in 1:4
        output[:,i,1] = out[i].B
        output[:,i,2] = out[i].mu_auto
        output[:,i,3] = out[i].mu_het
        output[:,i,4] = out[i].CD_auto
        output[:,i,5] = out[i].CD_tot
    end
    
    summ = zeros(nz, 3, 5)

    for i in 1:5
        summ[:,1,i] = mean(output[:,:,i], dims = 2)
        summ[:,2,i] = minimum(output[:,:,i], dims = 2)
        summ[:,3,i] = maximum(output[:,:,i], dims = 2)
    end
    
    return summ, output
end

nothing