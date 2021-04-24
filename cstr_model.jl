using JuMP
using Ipopt

println("tic");

# Setup
ntfe = 10
horizon = 15.0
tfe_width = horizon/ntfe
time = zeros(ntfe+1)
t0 = time[1]
for i = 0:ntfe
    time[i+1] = i*tfe_width
end

comp = ("A", "B")
stoich = Dict( [("A", -1), ("B", 1)] )

m = Model(Ipopt.Optimizer)

@variable(m, conc[time, comp])
@variable(m, dcdt[time, comp])

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
        dcdt[t, j] - stoich[j]*(conc[t, "A"] - conc[t, "B"]) == 0,
        )

fix(conc[t0, "A"], 1.0)
fix(conc[t0, "B"], 0.0)

optimize!(m)

println("toc");
