##### P2P Main File #####

##### activate environment #####
using Pkg
Pkg.activate(pwd())

##### packages #####
using JuMP

##### data read #####
data_file = "_input/input_data.xlsx";
@info "The data input file is set to" data_file
# data file
include("p2p_model-julia_data_read.jl");

##### MCP model #####

# model file
include("p2p_model-julia_mcp.jl");

# model run specific parameters
p_eex = 6.44;
p_dso = 3.48;
p_tso = 3.44;
p_eeg = 6.52;
p_td = 8.85;
@info "price parameters" p_eex p_dso p_tso p_eeg p_td

p_i = p_dso + p_tso + p_eeg + p_td;         # consumption from LBM
p_g = p_eex + p_dso + p_tso + p_eeg + p_td; # consumption from grid

## build and run P2P MCP model
 
setups = [1:6;];

for s in setups

    if s==1

        @info "1 - BAU Feed-in"

        p2p_mcp = build_p2p_mcp_model(N, T, p_g, p_i, p_d)

        fix.(p2p_mcp[:I], 0; force=true)      # no import
        fix.(p2p_mcp[:X], 0; force=true)      # no export
        fix.(p2p_mcp[:S_C], 0; force=true)    # no storage
        fix.(p2p_mcp[:S_D], 0; force=true)    # no storage
        fix.(p2p_mcp[:S], 0; force=true)      # no storage

    elseif s==2

        @info "2 - Local Sharing"
        p_i = p_dso + p_td; # consumption from LBM

        p2p_mcp = build_p2p_mcp_model(N, T, p_g, p_i, p_d)

        fix.(p2p_mcp[:F], 0; force=true)      # no feed-in
        fix.(p2p_mcp[:S_C], 0; force=true)    # no storage
        fix.(p2p_mcp[:S_D], 0; force=true)    # no storage
        fix.(p2p_mcp[:S], 0; force=true)      # no storage

    elseif s==3

        @info "3 - Home Storage"

        p2p_mcp = build_p2p_mcp_model(N, T, p_g, p_i, p_d)

        fix.(p2p_mcp[:I], 0; force=true)      # no import
        fix.(p2p_mcp[:X], 0; force=true)      # no export
        fix.(p2p_mcp[:F], 0; force=true)      # no feed-in

    elseif s==4

        @info "4 - Home Storage & Local Sharing"
        p_i = p_dso + p_td; # consumption from LBM

        p2p_mcp = build_p2p_mcp_model(N, T, p_g, p_i, p_d)

        fix.(p2p_mcp[:F], 0; force=true)      # no feed-in

    elseif s==5

        @info "5 - Current Regulatory Framework for 4"
        p_i = p_dso + p_tso + p_eeg + p_td; # consumption from LBM

        p2p_mcp = build_p2p_mcp_model(N, T, p_g, p_i, p_d)

        fix.(p2p_mcp[:F], 0; force=true)      # no feed-in

    elseif s==6

        @info "6 - Tech4All"
        p_i = p_dso + p_tso + 0.4 * p_eeg + p_td; # consumption from LBM
        p_d = p_d .+ 0.4 * p_eeg; # consumption from storage

        p2p_mcp = build_p2p_mcp_model(N, T, p_g, p_i, p_d)

        fix.(p2p_mcp[:F], 0; force=true)      # no feed-in

    end

    set_optimizer_attribute(p2p_mcp, "output", "no")
    set_optimizer_attribute(p2p_mcp, "output_options", "yes")
    set_optimizer_attribute(p2p_mcp, "output_errors", "yes")
    set_optimizer_attribute(p2p_mcp, "output_major_iterations", "yes")
    set_optimizer_attribute(p2p_mcp, "output_final_statistics", "yes")
    set_optimizer_attribute(p2p_mcp, "output_final_summary", "yes")
    set_optimizer_attribute(p2p_mcp, "cumulative_iteration_limit", 100000)
    set_optimizer_attribute(p2p_mcp, "minor_iteration_limit", 10000)
    set_optimizer_attribute(p2p_mcp, "convergence_tolerance", 1e-2)
    set_optimizer_attribute(p2p_mcp, "lemke_start", "always")

    optimize!(p2p_mcp)
    @show termination_status(p2p_mcp)
    export_parameters_mcp(p2p_mcp,s)

end