using JuMP
using Ipopt

println("Beginning CSTR model script");

function make_model(time, tfe_width)
    ntfe = length(time)-1
    comp = ["A", "B"]
    stoich = Dict( [("A", -1), ("B", 1)] )

    m = Model(Ipopt.Optimizer)

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
    	rate_gen[t, j] - stoich[j]*conc[t, "A"] == 0,
        )

    # Discretization equations
    @constraint(
            m,
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

    # Derivative equations
    @constraint(
            m,
            conc_diff_eqn[t=time, j=comp],
            dcdt[t, j] - (
                flow_in[t]*conc_in[t, j] -
                flow_out[t]*conc_out[t, j] +
                rate_gen[t, j]
                ) == 0,
            )

    return m
end

# Set-up time "set"
ntfe = 10
horizon = 60.0
tfe_width = horizon/ntfe
time = zeros(ntfe+1)
t0 = time[1]
for i = 0:ntfe
    time[i+1] = i*tfe_width
end

m = make_model(time, tfe_width)

for t=time
    fix(m[:conc_in][t, "A"], 5.0)
    fix(m[:conc_in][t, "B"], 0.01)
    if t < 10
        fix(m[:flow_in][t], 1.0)
    else
        fix(m[:flow_in][t], 5.0)
    end
end

fix(m[:conc][t0, "A"], 1.0)
fix(m[:conc][t0, "B"], 0.0)

optimize!(m)

println(value.(m[:conc]))

println("Finished CSTR model script");
