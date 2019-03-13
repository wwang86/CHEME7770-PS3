#ps3 part c
#Weiyao Wang

#given Flux.jl
include("Flux.jl")
using LinearAlgebra

#stoichiometric_matrix from part a)
Stoichiometric_matrix = [0 0 0 -1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
-1 0 0 1 -1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
-1 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 1 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0;
0 1 -1 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
-1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
0 0 1 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0;
0 0 0 1 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0;
1 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0 0;
1 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0 0 0;
0 0 0 0 -1 1 0 0 0 0 0 0 0 0 0 -1 0 0 0 0;
0 0 0 0 2 -2 0 0 0 0 0 0 0 0 0 0 -1 0 0 0;
0 0 0 0 1.5 -1.5 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
0 0 0 0 1.5 -1.5 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
0 0 0 0 -1.5 1.5 0 0 0 0 0 0 0 0 0 0 0 0 0 -1;
0 0 -1 0 -2 2 0 0 0 0 0 1 0 0 0 0 0 0 0 0]

#Define Kcat values
kcat1 = 203.0 #sec^-1
kcat2 = 34.5 #sec^-1
kcat3 = 249.0 #sec^-1
kcat4 = 88.1 #sec^-1
kcat5 = 13.7 #sec^-1
E = 0.01 #umol/gDW
water_fraction = 0.798 #Wet fraction of a cell
mass_cell = 2.3E-9 #g
vol_cell = 3.7E-13 #L
E1 = E / water_fraction * mass_cell / vol_cell  #umol/L

#Km values from Park et al. unit umol/L
Km_v5_arg = 4.4
Km_v5_NADPH = .3
Km_v4_orn = 360.0
Km_v3_arg = 1500.0
Km_v1_asp = 900.0
Km_v1_ATP = 51.0


#substrate concentration from Park et al. (others are assumed 0 if not found) unit umol/L
asp = 14900
orn = 4490.0
arg = 255.0
ATP = 4670.0
NADPH = 65.4

#Saturation
Saturation_1 = (asp/(asp+Km_v1_asp))*(ATP/(ATP+Km_v1_ATP))
Saturation_3 = (arg/(arg+Km_v3_arg))
Saturation_4 = (orn/(orn+Km_v4_orn))
Saturation_5 = (arg/(arg+Km_v5_arg))*(NADPH/(NADPH+Km_v5_NADPH))

#upper bound calculations
upper_bound_1 = kcat1 * Saturation_1 * E1
upper_bound_2 = kcat2 * E1
upper_bound_3 = kcat3 * Saturation_3 * E1
upper_bound_4 = kcat4 * Saturation_4 * E1
upper_bound_5 = kcat5 * Saturation_5* E1
upper_bound_b = 10000 / water_fraction * mass_cell / vol_cell /3600 #umol/L

Default_bounds = [0.0 upper_bound_1;
 0 upper_bound_2;
 0 upper_bound_3;
 0 upper_bound_4;
 0 upper_bound_5;
 0 upper_bound_5;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b;
 0 upper_bound_b]

Species_bounds_array=zeros(18,2);

Objective_coefficient = [0.0; 0; 0; 0; 0; 0; 0; 0; 0; -1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]

optimize = calculate_optimal_flux_distribution(Stoichiometric_matrix, Default_bounds_array, Species_bounds_array, Objective_coefficient_array)
r=3600*optimize[2];
urea_flux=r[10];
println("maximum urea flux=", urea_flux)
