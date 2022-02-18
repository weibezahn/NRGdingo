##### P2P MCP Model #####

using GAMS

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
function export_parameters_mcp(model,job,season,jid)
    # export parameter
    R_exp = DataFrame(value.(model[:R]), H);                        # renewables
    I_exp = DataFrame(value.(model[:I]), H);                        # import
    X_exp = DataFrame(value.(model[:X]), H);                        # export
    S_c_exp = DataFrame(value.(model[:S_c]), H);                    # charge
    S_d_exp = DataFrame(value.(model[:S_d]), H);                    # discharge
    S_exp = DataFrame(value.(model[:S]), H);                        # storage level
    P_P_exp = DataFrame(value.(model[:P_P]), H);
    P_LEM_exp = DataFrame(:value => [value.(model[:P_LEM])[x] for x in keys(value.(model[:P_LEM]))]);
    P_BAT_exp = DataFrame(value.(model[:P_BAT]), H);

    # Excel export
    XLSX.writetable(
        "_output/$jid-$job-s$season-results_mcp.xlsx",
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

## Solve
MOI.set(p2p_mcp, MOI.TimeLimitSec(), 100000)
optimize!(p2p_mcp)

## Solution status
@info termination_status(p2p_mcp)
@info solution_summary(p2p_mcp)
@info has_values(p2p_mcp)
@info solve_time(p2p_mcp)

## Post-processing
@info "exporting MCP results"
P_LEM_exp = export_parameters_mcp(p2p_mcp,job,season,jid)