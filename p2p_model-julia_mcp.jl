##### P2P MCP Model #####

using GAMS
# using PATHSolver

# PATHSolver.c_api_License_SetString("2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0")

# "Write Variable Output to CSV File"
# function write_output(job::String,variable::DataFrame,name::String)
#     if isdir("_output/$job") == false
#         mkdir("_output/$job")
#     end
#     if size(variable,2) == 1
#         CSV.write(string("_output/$job/$job-$name.csv"),variable,delim=";",header=[Symbol("N")])
#     else
#         CSV.write(string("_output/$job/$job-$name.csv"),variable,delim=";",header=[Symbol("H$i") for i in 1:size(variable,2)])
#     end
# end

## Functions

"""
Function to build the MCP model
"""
function build_p2p_mcp_model(T, N, p_i)

    ## Model
    p2p_mcp = Model();

    ## Variables
    @variables p2p_mcp begin
    R[t in T, n in N] >= 0
    # G[t in T, n in N] >= 0
    I[t in T, n in N] >= 0
    X[t in T, n in N] >= 0
    S_c[t in T, n in N] >= 0
    S_d[t in T, n in N] >= 0
    S[t in T, n in N] >= 0

    P_P[t in T, n in N]
    P_LEM[t in T]
    P_BAT[t in T, n in N]
    λ_c[t in T, n in N] >= 0
    λ_slb[t in T, n in N] >= 0
    λ_sub[t in T, n in N] >= 0
    λ_sc[t in T, n in N] >= 0
    λ_sd[t in T, n in N] >= 0
    end

    ## Constraints

    # ⟂ R
    @constraint(p2p_mcp, kkt1[t in T, n in N],
        p_lcoe[n] - P_P[t,n] + λ_c[t,n]  ⟂ R[t,n]
    )

    # ⟂ G
    # @constraint(p2p_mcp, kkt2[t in T, n in N],
    #     - P_P[t,n] + p_g ⟂ G[t,n]
    # )

    # ⟂ I
    @constraint(p2p_mcp, kkt3[t in T, n in N],
        p_i + P_LEM[t] - P_P[t,n] ⟂ I[t,n]
    )

    # ⟂ X
    @constraint(p2p_mcp, kkt4[t in T, n in N],
        - P_LEM[t] + P_P[t,n] ⟂ X[t,n]
    )

    # ⟂ S_c
    @constraint(p2p_mcp, kkt5[t in T, n in N],
        P_P[t,n] + η[n] * P_BAT[t,n] + λ_sc[t,n] ⟂ S_c[t,n]
    )

    # ⟂ S_d
    @constraint(p2p_mcp, kkt6[t in T, n in N],
        p_d[n] - P_P[t,n] - P_BAT[t,n] + λ_sd[t,n] ⟂ S_d[t,n]
    )

    # ⟂ S
    @constraint(p2p_mcp, kkt7a[t in T[1:end-1], n in N],
        - P_BAT[t,n] + P_BAT[t+1, n] - λ_slb[t,n] + λ_sub[t,n] ⟂ S[t,n]
    )
    @constraint(p2p_mcp, kkt7b[t in T[end], n in N],
        - P_BAT[t,n] - λ_slb[t,n] + λ_sub[t,n] ⟂ S[t,n]
    )

    # ⟂ λ_c
    @constraint(p2p_mcp, con1[t in T, n in N],
        res[t,n] - R[t,n] ⟂ λ_c[t,n]
    )

    # , P_P
    @constraint(p2p_mcp, con2[t in T, n in N],
        R[t,n] + I[t,n] + S_d[t,n] - X[t,n] - S_c[t,n] - dem[t,n] ⟂ P_P[t,n]
    )

    # ⟂ P_LEM
    @constraint(p2p_mcp, con3[t in T],
        sum(X[t,n] for n in N) - sum(ϕ * I[t,n] for n in N) ⟂ P_LEM[t]
    )

    # ⟂ λ_sd
    @constraint(p2p_mcp, con4[t in T, n in N],
        β[n] - S_d[t,n] ⟂ λ_sd[t,n]
    )

    # ⟂ λ_sc
    @constraint(p2p_mcp, con5[t in T, n in N],
        α[n] - S_c[t,n] ⟂ λ_sc[t,n]
    )

    # ⟂ λ_slb
    @constraint(p2p_mcp, con6[t in T, n in N],
        S[t,n] - s_min[n] ⟂ λ_slb[t,n]
    )

    # ⟂ λ_sub
    @constraint(p2p_mcp, con7[t in T, n in N],
        s_max[n] - S[t,n] ⟂ λ_sub[t,n]
    )

    # , P_BAT
    @constraint(p2p_mcp, con8a[t in T[1], n in N],
        S[t,n] - s_init - η[n] * S_c[t,n] + S_d[t,n] ⟂ P_BAT[t,n]
    )
    @constraint(p2p_mcp, con8b[t in T[2:end], n in N],
        S[t,n] - S[t-1,n] - η[n] * S_c[t,n] + S_d[t,n] ⟂ P_BAT[t,n]
    )

    # return model
    return p2p_mcp

end

"""
Function to export the results of the MCP model
"""
function export_parameters_mcp(model,job_id)
    # export parameter

    R_exp = DataFrame(value.(model[:R]), H);                        # renewables
    # write_output(job,out_R,"R")
    # G_exp = DataFrame(value.(G), H);
    # write_output(job,out_G,"G")
    I_exp = DataFrame(value.(model[:I]), H);                        # import
    # write_output(job,out_I,"I")
    X_exp = DataFrame(value.(model[:X]), H);                        # export
    # write_output(job,out_X,"X")
    S_c_exp = DataFrame(value.(model[:S_c]), H);                    # charge
    # write_output(job,out_S_c,"S_c")
    S_d_exp = DataFrame(value.(model[:S_d]), H);                    # discharge
    # write_output(job,out_S_d,"S_d")
    S_exp = DataFrame(value.(model[:S]), H);                        # storage level
    # write_output(job,out_S,"S")
    P_P_exp = DataFrame(value.(model[:P_P]), H);
    # write_output(job,out_P_P,"P_P")
    P_LEM_exp = DataFrame(:value => [value.(model[:P_LEM])[x] for x in keys(value.(model[:P_LEM]))]);
    # write_output(job,out_P_LEM,"P_LEM")
    P_BAT_exp = DataFrame(value.(model[:P_BAT]), H);
    # write_output(job,out_P_BAT,"P_BAT")

    # Excel export
    XLSX.writetable(
        "_output/results_mcp-$job.xlsx",
        Renewables = ( collect(DataFrames.eachcol(R_exp)), DataFrames.names(R_exp) ),
        Import = ( collect(DataFrames.eachcol(I_exp)), DataFrames.names(I_exp) ),
        Export = ( collect(DataFrames.eachcol(X_exp)), DataFrames.names(X_exp) ),
        Charge = ( collect(DataFrames.eachcol(S_c_exp)), DataFrames.names(S_c_exp) ),
        Discharge = ( collect(DataFrames.eachcol(S_d_exp)), DataFrames.names(S_d_exp) ),
        Storage_level = ( collect(DataFrames.eachcol(S_exp)), DataFrames.names(S_exp) ),
        P_P = ( collect(DataFrames.eachcol(P_P_exp)), DataFrames.names(P_P_exp) ),
        P_LEM = ( collect(DataFrames.eachcol(P_LEM_exp)), DataFrames.names(P_LEM_exp) ),
        P_BAT = ( collect(DataFrames.eachcol(P_BAT_exp)), DataFrames.names(P_BAT_exp) )
    )

    return P_LEM_exp

end

## Build P2P MCP model
@info "building MCP model with" p_i
p2p_mcp = build_p2p_mcp_model(T, N, p_i)

## GAMS

set_optimizer(p2p_mcp, GAMS.Optimizer)

set_optimizer_attribute(p2p_mcp, GAMS.License(), GAMS_lic_path)
set_optimizer_attribute(p2p_mcp, GAMS.ModelType(), "MCP")

## PATHSolver
# set_optimizer(p2p_mcp, PATHSolver.Optimizer)

# MOI.set(p2p_mcp, MOI.RawParameter("output"), "no")
# set_optimizer_attribute(p2p_mcp, "max_iter", 1000000)
# set_optimizer_attribute(p2p_mcp, "major_iteration_limit", 1000000)
# set_optimizer_attribute(p2p_mcp, "minor_iteration_limit", 1000000)
# set_optimizer_attribute(p2p_mcp, "cumulative_iteration_limit", 1000000)
# set_optimizer_attribute(p2p_mcp, "time_limit", 1000000)
# set_optimizer_attribute(p2p_mcp, "factorization_method", "blu_lusol")
# set_optimizer_attribute(p2p_mcp, "factorization_library_name", PATHSolver.LUSOL_LIBRARY_PATH)

## Solve
MOI.set(p2p_mcp, MOI.TimeLimitSec(), 100000)
optimize!(p2p_mcp)

## Dirty fix for GAMS model type
set_optimizer_attribute(p2p_mcp, GAMS.ModelType(), "MCP")
@info "solving" p2p_mcp
optimize!(p2p_mcp)

## Solution status
@show termination_status(p2p_mcp)
solution_summary(p2p_mcp)
termination_status(p2p_mcp)
has_values(p2p_mcp)
solve_time(p2p_mcp)

## Post-processing
@info "exporting MCP results"
P_LEM_exp = export_parameters_mcp(p2p_mcp,job)