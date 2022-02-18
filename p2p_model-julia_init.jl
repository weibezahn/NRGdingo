##### Julia P2P Model #####
# main file
# to be started as cluster job with tasks for seasons

##### activate environment #####
using Pkg
Pkg.activate(pwd())

##### read script input parameter #####
# (comment out for local use)
par_job = ARGS[1] # job name
par_jid = ENV["JOB_ID"] # job id
par_season = ENV["SGE_TASK_ID"] # task id (season)

##### packages #####
using JuMP

##### options #####
GAMS_lic_path = "gamslice.txt"; # specify path to GAMS license file
## local (comment out for cluster use) ##
#= job = "2021_12_10-pirgacha";
season = 1;
jid = 0; =#
## from script parameter (comment out for local use) ##
job = par_job;
jid = par_jid;
season = par_season;

##### data read #####
data_file = "_input/input_data-season_$season.xlsx";
@info "The data input file is set to" data_file
include("p2p_model-julia_data_read.jl")

##### MCP model #####
# model run specific parameters
p_i = 482.39;
@info "The MCP parameters are set to" p_i
# model file
include("p2p_model-julia_mcp.jl")

##### LP models #####
# model run specific parameters
p_cppa = 13.61;     # CPPA
p_i_I = 482.39;     # case I
p_i_II = 208.73;    # case II
p_i_IIa = 235.80;   # case IIa
@info "The LP parameters are set to" p_cppa p_i_I p_i_II p_i_IIa
# options
plem_load = false;   # true = load p_lem from file
# model file
include("p2p_model-julia_lp.jl")