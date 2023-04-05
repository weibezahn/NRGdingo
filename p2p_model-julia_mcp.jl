##### P2P MCP Model #####

using PATHSolver
PATHSolver.c_api_License_SetString("2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0");

"""
Function to build the MCP model
"""
function build_p2p_mcp_model(N, T, p_g, p_i, p_d)

    @info "building MCP model with" p_g p_i p_d

    ## model
    p2p_mcp = Model(PATHSolver.Optimizer)

    ## variables
    @variables p2p_mcp begin
        R[n in N, t in T] >= 0, (container=Array)
        G[n in N, t in T] >= 0, (container=Array)
        I[n in N, t in T] >= 0, (container=Array)
        X[n in N, t in T] >= 0, (container=Array)
        F[n in N, t in T] >= 0, (container=Array)
        S_C[n in N, t in T] >= 0, (container=Array)
        S_D[n in N, t in T] >= 0, (container=Array)
        S[n in N, t in T] >= 0, (container=Array)
     
        P_N[n in N, t in T], (container=Array)
        P_S[n in N, t in T], (container=Array)
        P_LBM[t in T], (container=Array)
        λ_res[n in N, t in T] >= 0, (container=Array)
        λ_slb[n in N, t in T] >= 0, (container=Array)
        λ_sub[n in N, t in T] >= 0, (container=Array)
        λ_α[n in N, t in T] >= 0, (container=Array)
        λ_β[n in N, t in T] >= 0, (container=Array)
     end;
     
     ## constraints

     # ⟂ R (Eq. A1)
     @constraint(p2p_mcp, kkt1[n in N, t in T],
         p_mc[n] - P_N[n,t] + λ_res[n,t] ⟂ R[n,t]
     );
     
     # ⟂ G (Eq. A2)
     @constraint(p2p_mcp, kkt2[n in N, t in T],
         p_g - P_N[n,t] ⟂ G[n,t]
     );
     
     # ⟂ I (Eq. A3)
     @constraint(p2p_mcp, kkt3[n in N, t in T],
         P_LBM[t] + p_i - P_N[n,t] ⟂ I[n,t]
     );
     
     # ⟂ X (Eq. A4)
     @constraint(p2p_mcp, kkt4[n in N, t in T],
         p_mc[n] - P_LBM[t] + λ_res[n,t] ⟂ X[n,t]
     );
     
     # ⟂ F (Eq. A5)
     @constraint(p2p_mcp, kkt5[n in N, t in T],
         p_mc[n] - p_fit[n] + λ_res[n,t] ⟂ F[n,t]
     );
     
     # ⟂ S_C (Eq. A6)
     @constraint(p2p_mcp, kkt6[n in N, t in T],
         p_mc[n] + λ_res[n,t] + η[n] * P_S[n,t] + λ_α[n,t] ⟂ S_C[n,t]
     );
     
     # ⟂ S_D (Eq. A7)
     @constraint(p2p_mcp, kkt7[n in N, t in T],
         p_d[n] - P_N[n,t] - P_S[n,t] + λ_β[n,t] ⟂ S_D[n,t]
     );
     
     # ⟂ S (Eq. A8)
     @constraint(p2p_mcp, kkt8a[n in N, t in T[1:end-1]],
         - P_S[n,t] + P_S[n,t+1] - λ_slb[n,t] + λ_sub[n,t] ⟂ S[n,t]
     );
     @constraint(p2p_mcp, kkt8b[n in N, t in T[end]],
         - P_S[n,t] - λ_slb[n,t] + λ_sub[n,t] ⟂ S[n,t]
     );
     
     # , P_N (Eq. A9)
     @constraint(p2p_mcp, con1[n in N, t in T],
         dem[n,t] - R[n,t] - G[n,t] - I[n,t] - S_D[n,t] ⟂ P_N[n,t]
     );
     
     # ⟂ λ_res (Eq. A10)
     @constraint(p2p_mcp, con2[n in N, t in T],
         res[n,t] - R[n,t] - X[n,t] - F[n,t] - S_C[n,t] ⟂ λ_res[n,t]
     );
     
     # , P_S (Eq. A11)
     @constraint(p2p_mcp, con3a[n in N, t in T[1]],
         s_init - S[n,t] + η[n] * S_C[n,t] - S_D[n,t] ⟂ P_S[n,t]
     );
     @constraint(p2p_mcp, con3b[n in N, t in T[2:end]],
         S[n,t-1] - S[n,t] + η[n] * S_C[n,t] - S_D[n,t] ⟂ P_S[n,t]
     );
     
     # ⟂ λ_slb (Eq. A12)
     @constraint(p2p_mcp, con4[n in N, t in T],
         S[n,t] - s_min[n] ⟂ λ_slb[n,t]
     );
     
     # ⟂ λ_sub (Eq. A13)
     @constraint(p2p_mcp, con5[n in N, t in T],
         s_max[n] - S[n,t] ⟂ λ_sub[n,t]
     );
     
     # ⟂ λ_α (Eq. A14)
     @constraint(p2p_mcp, con6[n in N, t in T],
         α[n] - S_C[n,t] ⟂ λ_α[n,t]
     );
     
     # ⟂ λ_β (Eq. A15)
     @constraint(p2p_mcp, con7[n in N, t in T],
         β[n] - S_D[n,t] ⟂ λ_β[n,t]
     );
     
     # , P_LBM (Eq. A23)
     @constraint(p2p_mcp, con8[t in T],
         sum(X[n,t] for n in N) - sum(I[n,t] for n in N) ⟂ P_LBM[t]
     );

    ## return model
    return p2p_mcp

end

## data output

using CSV

"""
Function to export the results of the MCP model
"""
function export_parameters_mcp(model,job)

    # export parameter
    R_exp = DataFrame(value.(model[:R])', H);                        # renewables
    G_exp = DataFrame(value.(model[:G])', H);                        # grid
    I_exp = DataFrame(value.(model[:I])', H);                        # import
    X_exp = DataFrame(value.(model[:X])', H);                        # export
    F_exp = DataFrame(value.(model[:F])', H);                        # feed-in
    S_C_exp = DataFrame(value.(model[:S_C])', H);                    # charge
    S_D_exp = DataFrame(value.(model[:S_D])', H);                    # discharge
    S_exp = DataFrame(value.(model[:S])', H);                        # storage level
    P_N_exp = DataFrame(value.(model[:P_N])', H);
    P_LBM_exp = DataFrame(:value => [value.(model[:P_LBM])[x] for x in keys(value.(model[:P_LBM]))]);
    P_S_exp = DataFrame(value.(model[:P_S])', H);

        # Excel export
        XLSX.writetable(
            "_output/$job-results_mcp.xlsx",
            Renewables = ( collect(DataFrames.eachcol(R_exp)), DataFrames.names(R_exp) ),
            Grid = ( collect(DataFrames.eachcol(G_exp)), DataFrames.names(G_exp) ),
            Import = ( collect(DataFrames.eachcol(I_exp)), DataFrames.names(I_exp) ),
            Export = ( collect(DataFrames.eachcol(X_exp)), DataFrames.names(X_exp) ),
            Feedin = ( collect(DataFrames.eachcol(F_exp)), DataFrames.names(F_exp) ),
            Charge = ( collect(DataFrames.eachcol(S_C_exp)), DataFrames.names(S_C_exp) ),
            Discharge = ( collect(DataFrames.eachcol(S_D_exp)), DataFrames.names(S_D_exp) ),
            Storage_level = ( collect(DataFrames.eachcol(S_exp)), DataFrames.names(S_exp) ),
            P_N = ( collect(DataFrames.eachcol(P_N_exp)), DataFrames.names(P_N_exp) ),
            P_LBM = ( collect(DataFrames.eachcol(P_LBM_exp)), DataFrames.names(P_LBM_exp) ),
            P_S = ( collect(DataFrames.eachcol(P_S_exp)), DataFrames.names(P_S_exp) )
        )
    
end