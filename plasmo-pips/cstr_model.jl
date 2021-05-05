#using JuMP
using Ipopt
using Plasmo
using MPIClusterManagers
using Distributed

@everywhere using Pkg
@everywhere Pkg.activate((@__DIR__))

@everywhere using Plasmo
@everywhere using PipsSolver

println("Beginning CSTR model script using Plasmo.jl");

function make_model(time, tfe_width, comp, stoich, k_rxn)
    nt = length(time)
    ntfe = nt-1

    graph = OptiGraph()
    @optinode(graph, nodes[time])

    for node in nodes
        @variable(node, conc[comp])
        @variable(node, dcdt[comp])

        @variable(node, flow_in)
        @variable(node, flow_out)
        @constraint(
        	node,
        	flow_eqn,
        	flow_in - flow_out == 0,
        	)

        @variable(node, conc_in[comp])
        @variable(node, conc_out[comp])
        @constraint(
        	node,
        	conc_out_eqn[j=comp],
        	conc_out[j] - conc[j] == 0,
        	)

        @variable(node, rate_gen[comp])
        @constraint(
        	node,
        	rate_eqn[j=comp],
        	rate_gen[j] - stoich[j]*k_rxn*conc["A"] == 0,
            )

        # Differential equations
        @constraint(
                node,
                conc_diff_eqn[j=comp],
                dcdt[j] - (
                    flow_in*conc_in[j] -
                    flow_out*conc_out[j] +
                    rate_gen[j]
                    ) == 0,
                )
    end

    # Discretization equations
    # 
    # These are linking equations, but nothing seems to change if
    # we delcare them with simply @constraint
    @linkconstraint(
            graph,
            dcdt_disc_eqn[t=time[2:ntfe+1], j=comp],
            nodes[t][:dcdt][j] -
            (
             nodes[t][:conc][j] -
             nodes[t-tfe_width][:conc][j]
            )/tfe_width == 0,
            # This implementation is unstable as it requires
            # (t-tfe_width) to give the correct answer. A more
            # stable implementation would find t-tfe_width with
            # a "binary-search-within-tolerance" function.
            #
            # Should probably define this constraint expression
            # with a function...
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
nodes = graph[:nodes]

println(graph)

conc_in = Dict( [("A", 5.0), ("B", 0.01)] )

for t=time
    fix(nodes[t][:conc_in]["A"], conc_in["A"])
    fix(nodes[t][:conc_in]["B"], conc_in["B"])
    if t == 0
        fix(nodes[t][:flow_in], 0.1)
    else
        fix(nodes[t][:flow_in], 1.0)
    end
end

fix(nodes[t0][:conc]["A"], 1.0)
fix(nodes[t0][:conc]["B"], 0.0)

# MPI code
manager = MPIManager(np=4)
addprocs(manager)

julia_workers = sort(collect(values(manager.mpi2j)))
remote_references = PipsSolver.distribute(graph, julia_workers, remote_name=:pipsgraph)

#ipopt = Ipopt.Optimizer
#optimize!(graph, ipopt)
#
#for t=time
#    println(nodevalue.(nodes[t][:conc]))
#    # Would like a syntax like:
#    # println(nodevalue.(nodes[:][:conc]))
#    #
#    # ... something like [:conc] "broadcasted"
#    # across all the nodes...
#end
#
#println("Finished CSTR model script");
