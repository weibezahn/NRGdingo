##### Julia P2P Model #####

## Activate environment
using Pkg
Pkg.activate(pwd())

## Packages
using JuMP

## Data read
data_file = "_input/input_data-mcp.xlsx";
@info "The data input file is set to" data_file
include("p2p_model-julia_data_read.jl")

## MCP Model
# model run specific parameters
p_i = 482.39;
#options
job = "Test_GAMS";
GAMS_lic_path = "/net/work/weibeza/gamslice---tub_generic---2021---exp_2022_04.txt";
include("p2p_model-julia_mcp.jl")

## LP Models
# model run specific parameters
p_cppa = 13.61;     # CPPA
p_i_II = 208.73;    # case II
p_i_IIa = 235.80;   # case IIa
@info "The LP parameters are set to" p_cppa p_i_II p_i_IIa
# options
plem_load = false;   # true = load p_lem from file
include("p2p_model-julia_lp.jl")