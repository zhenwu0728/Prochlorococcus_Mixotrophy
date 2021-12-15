# ##################################################################
# Julia code to model macromolecular composition, elemental composition
# and growth rate of photo-autotroph following Inomura et al (2020)
# Inomura, K., et al, "A mechanistic model of macromolecular ... ",
# Frontiers in Microbiology, 11, Article 86 (2020)
#
# Reads in physiological data from Liefer et al (2019) for comparison
# Liefer JD, Garg A, Fyfe MH, Irwin AJ, Benner I, Brown CM, Follows MJ, 
# Omta AW and Finkel ZV (2019) The Macromolecular Basis of Phytoplankton 
# C:N:P Under Nitrogen Starvation. 
# Front. Microbiol. 10:763. doi: 10.3389/fmicb.2019.00763
######################################################################
# 
#
# Mick Follows July/Aug 2020
# ##################################################################

# Packages
using Printf             # formatted printing
using Parameters         # pass many parameters between functions
#using Plots              # making plots, default backend GR
using CSV, DataFrames    # read CSV files into Dataframes
using Polynomials	 # fit and plot polynomial to data set
#using LaTeXStrings      # LaTeX typesetting in plots

# SET PARAMETERS for Inomura et al model
# use parameters.jl; @with_kw allows decoration with default values
@with_kw struct Phyto_params
# Params optimized for S. linearis: Inomura et al (2020), Table S5
	m::Float64 = 3.93e-1                # day-1
	nu_I_max::Float64 = 2.77e2          # mol C (mol C in Chl)-1 day-1
	A_I::Float64 = 8.63e-3              # micromol-1 m2 s
	A_Pho::Float64 = 1.60e1             # mol C (mol C in Chl)-1
	A_Bio::Float64 = 2.71e-1           # mol C (mol C)-1 day-1
	Q_C_Pro_Other::Float64 = 2.40e-1    # mol C (mol C)-1
	A_RNA_P::Float64 = 4.23e-3          # mol P (mol C)-1 day
	# NOTE A_Pho_PChl instead of Y_Pho_PChl (typo in Table S5)
	A_Pho_PChl::Float64 = 2.83e-2       # mol P (mol C in Chl)-1
	Q_P_Other0::Float64 = 6.53e-4       # mol P (mol C)-1
	Q_N_Sto_max::Float64 = 3.50e-2      # mol N (mol C)-1
	Q_C_Other0::Float64 = 1.82e-2       # mol C (mol C)-1
# Prescribed parameters: Inomura et al (2020) Supplementary Table 4
	E_Inomura::Float64 = 7.74e-1        # dimensionless
	Q_C_DNA::Float64 = 9.41e-4          # mol C (mol C)-1
	Q_Pmin_RNA::Float64 = 2.23e-4       # mol P (mol C)-1
	Q_P_max::Float64 = 2.23e-4          # mol P (mol C)-1
# Elemental ratios of macro-molecules: Inomura et al (2020) Table 1
	Y_Chl_CN::Float64 = 55.0/4.0
	Y_Chl_NC::Float64 = 1.0/Y_Chl_CN
	Y_Pro_CN::Float64 = 4.49/1.0
	Y_Pro_NC::Float64 = 1.0/Y_Pro_CN
	Y_RNA_CN::Float64 = 10.7/3.8
	Y_RNA_NC::Float64 = 1.0/Y_RNA_CN
	Y_RNA_CP::Float64 = 10.7/1.0
	Y_RNA_PC::Float64 = 1.0/Y_RNA_CP
	Y_RNA_NP::Float64 = 3.8/1.0
	Y_RNA_PN::Float64 = 1.0/Y_RNA_NP
	Y_DNA_CN::Float64 = 11.1/3.8
	Y_DNA_NC::Float64 = 1.0/Y_DNA_CN
	Y_DNA_CP::Float64 = 11.1/1.0
	Y_DNA_PC::Float64 = 1.0/Y_DNA_CP
	Y_DNA_NP::Float64 = 3.8/1.0
	Y_DNA_PN::Float64 = 1.0/Y_DNA_NP
	Y_Plip_CP::Float64 = 40.0/1.0
	Y_Plip_PC::Float64 = 1.0/Y_Plip_CP
	Y_Nstore_CN::Float64 = 2.0/1.0
	Y_Nstore_NC::Float64 = 1.0/Y_Nstore_CN
# some default parameter values (set mu and I in main program)
	E::Float64 = E_Inomura   # synth cost E = E_Inomura
	Q_N_Sto::Float64 = 0.0   # N storage, 0 by default
	Q_P_Sto::Float64 = 0.0   # P storage, 0 by default
end
pa = Phyto_params()
# END of Inomura et al parameters ......................

# provide "default" values for growth rate and light
mu = 0.4       # growth rate,   day-1
I = 144.0      # light intensity,   micromol m-2 s-1
A_Chl = 0.0    # default
B_Chl = 0.0    # default

# Functions: solutions from Inomura et al (2020) Supp. Table 3

# nu_I(I), photosynthesis vs light. Evaluate A_Chl and B_Chl........
function photosynth(I,pa) #.........................................
 	@unpack nu_I_max, A_I, E, m = pa   # get params used from pa
	nu_I = nu_I_max*(1.0 - exp(-A_I*I))
	A_Chl = (1.0 + E)/nu_I
	B_Chl = m / nu_I
	return nu_I, A_Chl, B_Chl
end #...............................................................

# mu_max_I, max light-limited growth rate for given I (day-1) ......
function mu_max_I_eval(A_Chl,B_Chl,pa) #............................
  	@unpack Y_RNA_CP, A_RNA_P, A_Pho, A_Bio, Y_Plip_CP, A_Pho_PChl = pa 
  	@unpack Q_C_Pro_Other, Q_C_Other0, Q_Pmin_RNA, Q_C_DNA = pa 
	# coefficients of quadratic
	a_M = Y_RNA_CP*A_RNA_P*(A_Pho*A_Chl + A_Bio)
 	b_M = (  (1.0 + A_Pho + Y_Plip_CP*A_Pho_PChl)*A_Chl + A_Bio
                  + Y_RNA_CP*A_RNA_P*(A_Pho*B_Chl + Q_C_Pro_Other)  )
	Q_C_Other = Q_C_Other0 + Q_C_DNA + Q_C_Pro_Other
 	c_M = ( (1.0 + A_Pho + Y_Plip_CP*A_Pho_PChl)*B_Chl
                  + Q_C_Other + Y_RNA_CP*Q_Pmin_RNA - 1.0     )
	# evaluate root of quadratic
	mu_max_I = (- b_M + (b_M*b_M - 4.0*a_M*c_M)^0.5 )/(2.0*a_M)
	return mu_max_I
end #...............................................................

# Chlorophyll to carbon ratio (moles C in Chl (mol C)-1) ...........
function ChlC(mu,A_Chl,B_Chl,pa) #..................................
#  	@unpack mu = pa
	Q_C_Chl =  A_Chl*mu + B_Chl
	return Q_C_Chl
end #...............................................................

# NC = N:C - nitrogen:carbon ratio .................................
function NC_eval(mu,A_Chl,B_Chl,pa) #................................
 	@unpack Y_RNA_CP, A_RNA_P, A_Pho, A_Bio, Y_Plip_CP, A_Pho_PChl = pa 
 	@unpack Y_RNA_NP, Q_C_Pro_Other, Q_C_Other0, Q_Pmin_RNA = pa 
 	@unpack Y_Chl_NC, Y_Pro_NC, Y_Chl_NC, Y_DNA_NC, Q_C_DNA = pa 
 	@unpack Q_N_Sto = pa 
	# coefficients of polynomial
	a_N = Y_RNA_NP*A_RNA_P*(A_Pho*A_Chl + A_Bio)
	b_N = (  (Y_Chl_NC*A_Chl + Y_Pro_NC*(A_Bio + A_Pho*A_Chl))
       	        + Y_RNA_NP*A_RNA_P*(A_Pho*B_Chl + Q_C_Pro_Other) )
	c_N = ( Y_Chl_NC*B_Chl + Y_Pro_NC*(A_Pho*B_Chl + Q_C_Pro_Other)
                + Y_RNA_NP*Q_Pmin_RNA + Y_DNA_NC*Q_C_DNA + Q_N_Sto )
	#println("c_N = ",c_N)
	# evaluate polynomial
#	NC(mu) = a_N*mu*mu + b_N*mu + c_N
	NC = a_N*mu*mu + b_N*mu + c_N
	return NC   
end #...............................................................

# PC = P:C - phosphorus:carbon ratio ...............................
function PC_eval(mu,A_Chl,B_Chl,pa) #................................
 	@unpack A_RNA_P, A_Pho, A_Bio, Y_Plip_CP, A_Pho_PChl, Q_C_DNA = pa 
 	@unpack Y_DNA_PC, Q_C_Pro_Other, Q_P_Other0, Q_Pmin_RNA = pa 
	@unpack Q_P_Sto = pa 
	# coefficients of polynomial
	a_P = A_RNA_P*(A_Bio + A_Pho*A_Chl)
	b_P = A_RNA_P*(A_Pho*B_Chl + Q_C_Pro_Other) + A_Pho_PChl*A_Chl
	c_P = ( Q_Pmin_RNA + Y_DNA_PC*Q_C_DNA + A_Pho_PChl*B_Chl
                 + Q_P_Other0 + Q_P_Sto )
	# evaluate polynomial
	# PC(mu) = a_P*mu*mu + b_P*mu + c_P
	PC = a_P*mu*mu + b_P*mu + c_P
	return PC
end #...............................................................

# Macro-molecular fractions of cellular carbon .....................
function MM_eval(mu,Q_C_Chl,pa) #...................................
 	@unpack A_RNA_P, A_Pho, A_Bio, Q_C_DNA, Q_C_Pro_Other = pa 
 	@unpack Y_Plip_CP, A_Pho_PChl, Y_RNA_CP, Q_Pmin_RNA = pa 
 	@unpack Q_C_Other0 = pa 
	# from protein, RNA, DNA, lipids, carbohydrate, Chl, Nstore
	#Q_C_DNA = Q_C_DNA      # prescribed param based on genome size
	#Q_C_Chl = Q_C_Chl     # solved using function ChlC 
	Q_C_Nstore = 0.0       # assuming no N-storage for now
	Q_C_Pro = A_Pho*Q_C_Chl + A_Bio*mu + Q_C_Pro_Other
	Q_C_RNA = (Q_Pmin_RNA + A_RNA_P*mu*Q_C_Pro)*Y_RNA_CP 
	Q_C_PlipThy = Q_C_Chl*A_Pho_PChl*Y_Plip_CP  
	# eval total C store: Q_C_store = Q_C_CarbStore +  Q_C_LipStore
	Q_C_store = 1.0 - (Q_C_Chl + Q_C_Pro + Q_C_RNA + Q_C_DNA + 
		 Q_C_PlipThy + Q_C_Nstore + Q_C_Other0) 
	# in order to evaluate carbs and lipid contributions separately
	# must define fraction of Q_C_Other0 and Q_C_store in each
	# - rather arbitrary at half and half right now...
	frac_carb = 0.5 
	Q_C_CarbOther = Q_C_Other0 * frac_carb
	Q_C_LipOther = Q_C_Other0 * (1.0 - frac_carb)
	Q_C_CarbStore = Q_C_store * frac_carb
	Q_C_LipStore = Q_C_store * (1.0 - frac_carb)
	Q_C_lip = Q_C_PlipThy + Q_C_LipStore + Q_C_LipOther
	Q_C_carb = Q_C_CarbStore + Q_C_CarbOther
	# Evaluate protein allocation contributions
	Q_C_Pro_Pho = A_Pho*Q_C_Chl
	Q_C_Pro_Bio = A_Bio*mu 
	Q_C_Pro_Other1 = Q_C_Pro_Other
	# Q_C_Chl and Q_C_DNA already evaluated
	# return other contributions
	#return Q_C_Pro, Q_C_RNA, Q_C_Nstore, Q_C_lip, Q_C_carb
	return Q_C_Pro, Q_C_RNA, Q_C_Nstore, Q_C_lip, Q_C_carb, Q_C_Pro_Pho, Q_C_Pro_Bio, Q_C_Pro_Other1
end #...............................................................