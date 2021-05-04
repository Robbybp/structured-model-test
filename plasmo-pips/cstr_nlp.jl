include("cstr_model.jl")

nodes = graph[:nodes]

for t=time
    flow_in = nodes[t][:flow_in]
    if t != t0
        unfix(flow_in)
        set_lower_bound(flow_in, 0.0)
    end
end

F_tgt = 1.0
conc_setpoint = Dict()
conc_setpoint["A"] = F_tgt*conc_in["A"]/(F_tgt - stoich["A"])
conc_setpoint["B"] = conc_in["B"] + stoich["B"]*conc_setpoint["A"]/F_tgt

sp_error = @expression(
        graph,
        sp_error[t=time, j=comp],
        nodes[t][:conc][j] .- conc_setpoint[j],
)

obj = @objective(graph, Min, sum(sp_error.^2))

ipopt = Ipopt.Optimizer
optimize!(graph, ipopt)

for t in time
    println(nodevalue.(nodes[t][:conc]))
    println(nodevalue.(nodes[t][:flow_in]))
end
