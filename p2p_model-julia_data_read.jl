##### P2P Data Read #####

using DataFrames
# using CSV
using XLSX

## Data read-in from XLSX file
demand_df = DataFrame(XLSX.readtable(data_file, "Demand")...);
lcoe_df = DataFrame(XLSX.readtable(data_file, "LCOE")...);
res_df = DataFrame(XLSX.readtable(data_file, "Renewables")...);
bat_df = DataFrame(XLSX.readtable(data_file, "Storage")...);

## Data read-in from CSV files
# demand_df = CSV.read("_input/data_demand.csv", DataFrame);
# lcoe_df = CSV.read("_input/data_lcoe.csv", DataFrame);
# res_df = CSV.read("_input/data_production.csv", DataFrame);
# bat_df = CSV.read("_input/data_battery.csv", DataFrame);

## Generation and demand data
res = Matrix(res_df[!,2:end]);      # renewable production of house n
dem = Matrix(demand_df[!,2:end]);   # demand of house n

## Cost data
p_lcoe = Matrix(lcoe_df[!,2:end]);      # LCOE of house n

# p_eex = 6.44
# p_tso = 3.44
# p_eeg = 6.52
# p_dso = 3.48
# p_td = 8.85

# p_g = p_eex + p_tso + p_eeg + p_dso + p_td
# p_i = p_tso + p_eeg + p_dso + p_td

p_d = bat_df[:,:p_d];       # price of discharging battery of house n

## Storage data
α = bat_df[:,:alpha];       # maximum charge rate of battery of house n
β = bat_df[:,:beta];        # maximum discharge rate of battery of house n
η = bat_df[:,:eta];         # battery discharging efficiency of house n
s_min = bat_df[:,:s_min];   # lower storage bound of house n
s_max = bat_df[:,:s_max];   # upper storage bound of house n
s_init = 0;                 # initial storage level

## Transmission data
ϕ = 0.93;                   # transmission loss

## Sets
T = [1:336...];             # hours t in time horizon T
N = [1:size(dem,2)...];     # houses n in community N
H = ["H$(N[i])" for i in 1:length(N)];

@info "The model parameter scope is set to" T N