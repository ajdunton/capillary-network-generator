# Author: Aaron Dunton (ajdunton@gmail.com)

# %% Manual Inputs
soil_HSG = "B" 

storm_intensity = 2.07 #in/hr

# %% Method Options
# Wastewater flow estimation: 'population' or 'buildings'
ww_method = 'population'

# Weights for spanning tree: 'hybrid-slope' or 'distance' or 'hybrid-elev'
weight_method = 'hybrid-elev'

# %% Parameter Definitions
# Parameters for determining nodes
max_drain_spacing = 400
maxratio = 1.5
maxdist = 150
max_node_spacing = 100

# Parameters for evaluating weight
s_mu = 0.002
s_sig = 0.004

# Hydrologic parameters
n_imp = 0.013
n_perv = 0.05
d_imp = 0.0625
d_perv = 0.25
if soil_HSG == "A": #Based on loamy sand, A.2 of SWMM User's Manual
    sat_hyd_cond = 1.18 
    suc_head = 2.4 
    porosity = 0.437 
    field_cap = 0.105 
    wilt_pt = 0.047
elif soil_HSG == 'B': #Based on loam, A.2 of SWMM User's Manual
    sat_hyd_cond = 0.13
    suc_head = 3.5
    porosity = 0.463
    field_cap = 0.232
    wilt_pt = 0.116
elif soil_HSG == 'C': #Based on sandy clay loam, A.2 of SWMM User's Manual
    sat_hyd_cond = 0.06
    suc_head = 8.66
    porosity = 0.398
    field_cap = 0.244
    wilt_pt = 0.136
elif soil_HSG == 'D': #Based on sandy clay, A.2 of SWMM User's Manual
    sat_hyd_cond = 0.02
    suc_head = 9.45
    porosity = 0.43
    field_cap = 0.321
    wilt_pt = 0.221

# Parameters for wastewater flow estimation
ww_comm = 2000/43560
ww_ind = 10000/43560
factor_ci = 2
factor_designpop = 2
ww_percap = 100
design_pop_single = 3.7
design_pop_multi = 3
ww_hotel = 130 

# Minimum pipe depth
min_depth = 5