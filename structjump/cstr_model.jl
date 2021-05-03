using JuMP
using StructJuMP
using Ipopt

println("Beginning structured CSTR model Julia script");

function make_model(time, tfe_width, comp, stoich, k_rxn)
    ntfe = length(time)-1

    time_index_map = Dict([
            (t, i) for (i, t) in enumerate(time)
            ])
    m = StructuredModel(num_scenarios=length(time))
    #m = Model()

    # It appears that any variable that appears in "linking constraints"
    # must be declared in the root model.
    @variable(m, conc[time, comp])
    @variable(m, dcdt[time, comp])
    for (i, t) in enumerate(time[2:ntfe+1])
        idx = i+1
        for j in comp
            # This raises a mysterious error:
            # type Int64 has no field value
            # Works with Model... is this a bug?
            @constraint(
                    m,
                    dcdt[t, j] -
                    (conc[t, j] - conc[t-tfe_width, j])/tfe_width == 0
                   )
        end
    end
    #@constraint(
    #        m,
    #        dcdt_disc_eqn[t=time[2:ntfe+1], j=comp],
    #        dcdt[t, j] -
    #        (conc[t, j] - conc[t-tfe_width, j])/tfe_width == 0,
    #        # This implementation is unstable as it requires
    #        # (t-tfe_width) to give the correct answer. A more
    #        # stable implementation would find t-tfe_width with
    #        # a "binary-search-within-tolerance" function.
    #        #
    #        # Should probably define this constraint expression
    #        # with a function...
    #        )

    for i in 1:ntfe
        t = time[i]
        bl = StructuredModel(parent=m, id=i)
        _conc = conc[t, :]
        _dcdt = dcdt[t, :]

        @variable(bl, flow_in)
        @variable(bl, flow_out)
        @constraint(
        	bl,
        	flow_eqn,
        	flow_in - flow_out == 0,
        	)

        @variable(bl, conc_in[comp])
        @variable(bl, conc_out[comp])
        @constraint(
        	bl,
        	conc_out_eqn[j=comp],
        	conc_out[j] - _conc[j] == 0,
        	)

        @variable(bl, rate_gen[comp])
        @constraint(
        	bl,
        	rate_eqn[j=comp],
        	rate_gen[j] - stoich[j]*k_rxn*_conc["A"] == 0,
            )

        # Differential equations
        @constraint(
                bl,
                conc_diff_eqn[j=comp],
                _dcdt[j] - (
                    flow_in*conc_in[j] -
                    flow_out*conc_out[j] +
                    rate_gen[j]
                    ) == 0,
                )
    end

    return m
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

m = make_model(time, tfe_width, comp, stoich, k_rxn)

println(m)

#conc_in = Dict( [("A", 5.0), ("B", 0.01)] )
#
#for t=time
#    fix(m[:conc_in][t, "A"], conc_in["A"])
#    fix(m[:conc_in][t, "B"], conc_in["B"])
#    if t == 0
#        fix(m[:flow_in][t], 0.1)
#    else
#        fix(m[:flow_in][t], 1.0)
#    end
#end
#
#fix(m[:conc][t0, "A"], 1.0)
#fix(m[:conc][t0, "B"], 0.0)

#optimize!(m)

#conc = m[:conc]
#println(value(conc))

println("Finished structured CSTR model Julia script");
