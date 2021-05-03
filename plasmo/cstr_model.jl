#using JuMP
using Ipopt
using Plasmo

println("Beginning CSTR model script using Plasmo.jl");

function make_model(time, tfe_width, comp, stoich, k_rxn)
    ntfe = length(time)-1

    #m = Model(Ipopt.Optimizer)
    graph = OptiGraph()
    @optinode(graph, m)

    @variable(m, conc[time, comp])
    @variable(m, dcdt[time, comp])

    @variable(m, flow_in[time])
    @variable(m, flow_out[time])
    @constraint(
    	m,
    	flow_eqn[t=time],
    	flow_in[t] - flow_out[t] == 0,
    	)

    @variable(m, conc_in[time, comp])
    @variable(m, conc_out[time, comp])
    @constraint(
    	m,
    	conc_out_eqn[t=time, j=comp],
    	conc_out[t, j] - conc[t, j] == 0,
    	)

    @variable(m, rate_gen[time, comp])
    @constraint(
    	m,
    	rate_eqn[t=time, j=comp],
    	rate_gen[t, j] - stoich[j]*k_rxn*conc[t, "A"] == 0,
        )

    # Discretization equations
    @constraint(
            graph,
            dcdt_disc_eqn[t=time[2:ntfe+1], j=comp],
            dcdt[t, j] -
                (conc[t, j] - conc[t-tfe_width, j])/tfe_width == 0,
            # This implementation is unstable as it requires
            # (t-tfe_width) to give the correct answer. A more
            # stable implementation would find t-tfe_width with
            # a "binary-search-within-tolerance" function.
            #
            # Should probably define this constraint expression
            # with a function...
            )

    # Differential equations
    @constraint(
            m,
            conc_diff_eqn[t=time, j=comp],
            dcdt[t, j] - (
                flow_in[t]*conc_in[t, j] -
                flow_out[t]*conc_out[t, j] +
                rate_gen[t, j]
                ) == 0,
            )

    return graph
end

# Set-up time "set"
ntfe = 10
horizon = 10.0
tfe_width = horizon/ntfe
time = zeros(ntfe+1)
t0 = time[1]
for i = 0:ntfe
    time[i+1] = i*tfe_width
end

comp = ["A", "B"]
stoich = Dict( [("A", -1), ("B", 1)] )
k_rxn = 1.0

graph = make_model(time, tfe_width, comp, stoich, k_rxn)
m = find_node(graph, 1)

println(graph)

conc_in = Dict( [("A", 5.0), ("B", 0.01)] )

for t=time
    fix(m[:conc_in][t, "A"], conc_in["A"])
    fix(m[:conc_in][t, "B"], conc_in["B"])
    if t == 0
        fix(m[:flow_in][t], 0.1)
    else
        fix(m[:flow_in][t], 1.0)
    end
end

fix(m[:conc][t0, "A"], 1.0)
fix(m[:conc][t0, "B"], 0.0)

ipopt = Ipopt.Optimizer
optimize!(graph, ipopt)

#println(value.(m[:conc]))

println("Finished CSTR model script");
