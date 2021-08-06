##### P2P LP Model #####

using Clp

## Functions

"""
Function to build the LP model with CPPA
"""
function build_p2p_cppa_lp_model(T, N, p_cppa, p_i, p_lem)
    
    ## Model
    p2p_lp = Model()

    ## Variables
    @variables p2p_lp begin
        R[t in T, n in N] >= 0      # renewables production of house n in time step t
        F[t in T, n in N] >= 0      # feed into the national grid of house n in time step t
        I[t in T, n in N] >= 0      # P2P electricity purchase of house n in time step t
        X[t in T, n in N] >= 0      # P2P electricity sale of house n in time step t
        S_c[t in T, n in N] >= 0    # charge of battery of house n in time step t
        S_d[t in T, n in N] >= 0    # discharge of battery of house n in time step t
        S[t in T, n in N] >= 0      # energy storage level in battery of house n in time step t
        # S_d[t in T[1],n in N] == 0)   # fix storage leven/discharge in first time step to 0?
    end

    ## Objective function
    @objective p2p_lp Min begin
        sum(p_lcoe[n] * R[t,n] + I[t,n] * (p_lem[t] + p_i) + S_d[t,n] * p_d[n] - X[t,n] * p_lem[t] - p_cppa * F[t,n] for t in T, n in N)
    end

    ## Constraints
    #  energy balance
    @constraint(p2p_lp, nrg_bal[t in T, n in N],
        R[t,n] + S_d[t,n] + I[t,n] == dem[t,n] + S_c[t,n] + X[t,n] + F[t,n]
    )

    # renewable production constraint
    @constraint(p2p_lp, r_con[t in T,n in N],
        R[t,n] <= res[t,n]
    )

    # node balance taking into account trade efficiency
    @constraint(p2p_lp, n_bal[t in T],
        sum(X[t,n] for n in N) == ϕ * sum(I[t,n] for n in N)
    )

    # community balance
    @constraint(p2p_lp, com_bal,
        sum(dem[t,n] + F[t,n] for t in T, n in N) <= sum(R[t,n] for t in T, n in N)
    )

    # energy storage level above lower bound
    @constraint(p2p_lp, s_min_con[t in T,n in N],
        S[t,n] >= s_min[n]
    )

    # energy storage level below upper bound
    @constraint(p2p_lp, s_max_con[t in T,n in N],
        S[t,n] <= s_max[n]
    )

    # charge rate below maximum charge rate
    @constraint(p2p_lp, alpha_con[t in T,n in N],
        S_c[t,n] <= α[n]
    )

    # discharge rate below maximum discharge rate
    @constraint(p2p_lp, beta_con[t in T,n in N],
        S_d[t,n] <= β[n]
    )

    # intertemporal energy storage level in battery
    @constraint(p2p_lp, s_con_1[t in T[1],n in N],
        S[t,n] == s_init + η[n] * S_c[t,n] - S_d[t,n]
    )
    @constraint(p2p_lp, s_con[t in T[2:end],n in N],
        S[t,n] == S[t-1,n] + η[n] * S_c[t,n] - S_d[t,n]
    )

    # return model
    return p2p_lp

end

"""
Function to export the results of the LP model
"""
function export_parameters_lp(model,p_i,job_id)
    
    # export parameter
    
    C_exp = [sum(p_lcoe[n] * value(model[:R][t,n]) + value(model[:I][t,n]) * (p_lem[t] + p_i) + value(model[:S_d][t,n]) * p_d[n] - value(model[:X][t,n]) * p_lem[t] - p_cppa * value(model[:F][t,n]) for t in T) / 100 for n in N]
    C_exp = DataFrame([H,C_exp], ["house","objective"])     # objective
    
    R_exp = DataFrame(value.(model[:R]), H)                 # renewables
    I_exp = DataFrame(value.(model[:I]), H)                 # import
    X_exp = DataFrame(value.(model[:X]), H)                 # export
    S_c_exp = DataFrame(value.(model[:S_c]), H)             # charge
    S_d_exp = DataFrame(value.(model[:S_d]), H)             # discharge
    S_exp = DataFrame(value.(model[:S]), H)                 # storage level

    # Excel export

    XLSX.writetable(
        "_output/results_lp-$job_id.xlsx",
        Objective = ( collect(DataFrames.eachcol(C_exp)), DataFrames.names(C_exp) ),
        Renewables = ( collect(DataFrames.eachcol(R_exp)), DataFrames.names(R_exp) ),
        Import = ( collect(DataFrames.eachcol(I_exp)), DataFrames.names(I_exp) ),
        Export = ( collect(DataFrames.eachcol(X_exp)), DataFrames.names(X_exp) ),
        Charge = ( collect(DataFrames.eachcol(S_c_exp)), DataFrames.names(S_c_exp) ),
        Discharge = ( collect(DataFrames.eachcol(S_d_exp)), DataFrames.names(S_d_exp) ),
        Storage_level = ( collect(DataFrames.eachcol(S_exp)), DataFrames.names(S_exp) )
    )

end

## Load p_lem from MCP model run
if plem_load == true
    ## load p_lem from file
    mcp_results_file = "_output/results_mcp-$job.xlsx"
    plem_df = DataFrame(XLSX.readtable(mcp_results_file, "P_LEM")...);
    # plem_df = CSV.read("_input/data_plem.csv", DataFrame);
    p_lem = Matrix(plem_df);
    @info "p_lem read from file" mcp_results_file
else
    ## load p_lem from cache
    p_lem = Matrix(P_LEM_exp)
    @info "p_lem read from variable" P_LEM_exp
end

## Build CPPA LP models
@info "building LP model case II with" p_cppa p_i_II
p2p_lp_II = build_p2p_cppa_lp_model(T, N, p_cppa, p_i_II, p_lem)
@info "building LP model case IIa with" p_cppa p_i_IIa
p2p_lp_IIa = build_p2p_cppa_lp_model(T, N, p_cppa, p_i_IIa, p_lem)

## Solve
set_optimizer(p2p_lp_II, Clp.Optimizer)
set_optimizer(p2p_lp_IIa, Clp.Optimizer)
@info "solving" p2p_lp_II
optimize!(p2p_lp_II)
@info "solving" p2p_lp_IIa
optimize!(p2p_lp_IIa)

## Post-processing
@info "exporting LP results"
export_parameters_lp(p2p_lp_II,p_i_II,"$job-II")
export_parameters_lp(p2p_lp_IIa,p_i_IIa,"$job-IIa")