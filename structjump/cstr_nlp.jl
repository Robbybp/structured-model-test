include("cstr_model.jl")

flow_in = m[:flow_in]
conc = m[:conc]

for t=time
    if t != t0
        unfix(flow_in[t])
        set_lower_bound(flow_in[t], 0.0)
    end
end

F_tgt = 1.0
conc_setpoint = Dict()
conc_setpoint["A"] = F_tgt*conc_in["A"]/(F_tgt - stoich["A"])
conc_setpoint["B"] = conc_in["B"] + stoich["B"]*conc_setpoint["A"]/F_tgt

sp_error = @expression(
        m,
        sp_error[t=time, j=comp],
        conc[t, j] .- conc_setpoint[j],
)

obj = @objective(m, Min, sum(sp_error.^2))

optimize!(m)

println(value.(conc))
println(value.(flow_in))
