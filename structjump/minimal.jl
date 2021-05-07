using JuMP
using StructJuMP

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
ntfe = length(time)-1

time_index_map = Dict([
        (t, i) for (i, t) in enumerate(time)
        ])
m = StructuredModel(num_scenarios=length(time))
#m = Model()

# It appears that any variable that appears in "linking constraints"
# must be declared in the root model.
@variable(m, conc[time])
@variable(m, dcdt[time])

_conc = Dict([(t, conc[t]) for t in time])
_dcdt = Dict([(t, dcdt[t]) for t in time])

@constraint(m,
            eq[t=time],
            _dcdt[t] - _conc[t] == 0,
           )

#@constraint(
#      m,
#      disc_eq[t=time[2:ntfe+1], j=comp],
#      dcdt[t, j] -
#      (conc[t, j] - conc[t-tfe_width, j])/tfe_width == 0,
#      )


println(m)
