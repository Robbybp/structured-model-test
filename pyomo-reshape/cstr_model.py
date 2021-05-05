import pyomo.environ as pyo
import pyomo.dae as dae
from pyomo.dae.flatten import flatten_dae_components


solver = pyo.SolverFactory("ipopt")


def _flow_eqn_rule(m, t):
    return m.flow_in[t] - m.flow_out[t] == 0


def _conc_out_eqn_rule(m, t, j):
    return m.conc[t, j] - m.conc_out[t, j] == 0


def _rate_eqn_rule(m, t, j):
    return m.rate_gen[t, j] - m.stoich[j]*m.k_rxn*m.conc[t, "A"] == 0


def _conc_diff_eqn_rule(m, t, j):
    return m.dcdt[t, j] - (
            m.flow_in[t]*m.conc_in[t, j] -
            m.flow_out[t]*m.conc_out[t, j] +
            m.rate_gen[t, j]
            ) == 0


def make_model(horizon=10.0):
    m = pyo.ConcreteModel()
    m.comp = pyo.Set(initialize=["A", "B"])
    m.time = dae.ContinuousSet(initialize=[0, horizon])
    time = m.time
    comp = m.comp

    m.stoich = pyo.Param(m.comp, initialize={"A": -1, "B": 1}, mutable=True)
    m.k_rxn = pyo.Param(initialize=1.0, mutable=True)

    m.conc = pyo.Var(m.time, m.comp)
    m.dcdt = dae.DerivativeVar(m.conc, wrt=m.time)

    m.flow_in = pyo.Var(time)
    m.flow_out = pyo.Var(time)
    m.flow_eqn = pyo.Constraint(time, rule=_flow_eqn_rule)

    m.conc_in = pyo.Var(time, comp)
    m.conc_out = pyo.Var(time, comp)
    m.conc_out_eqn = pyo.Constraint(time, comp, rule=_conc_out_eqn_rule)

    m.rate_gen = pyo.Var(time, comp)
    m.rate_eqn = pyo.Constraint(time, comp, rule=_rate_eqn_rule)

    m.conc_diff_eqn = pyo.Constraint(time, comp, rule=_conc_diff_eqn_rule)

    return m


def initialize_model(m, ntfe=10):
    disc = pyo.TransformationFactory("dae.finite_difference")
    disc.apply_to(m, wrt=m.time, nfe=ntfe, scheme="BACKWARD")

    t0 = m.time.first()

    m.conc_in[:, "A"].fix(5.0)
    m.conc_in[:, "B"].fix(0.01)

    m.flow_in[:].fix(1.0)
    m.flow_in[t0].fix(0.1)

    m.conc[t0, "A"].fix(1.0)
    m.conc[t0, "B"].fix(0.0)


def reshape(m):
    time = m.time

    # Create block to hold time-decomposition
    m.time_block = pyo.Block(time)

    # Identify time-indexed components
    scalar_vars, dae_vars = flatten_dae_components(m, time, pyo.Var)
    scalar_cons, dae_cons = flatten_dae_components(m, time, pyo.Constraint)

    for t in time:
        b = m.time_block[t]
        var_list = []
        b.vars = pyo.Reference([var[t] for var in dae_vars])
        con_list = []
        for con in dae_cons:
            try:
                condata = con[t]
                con_list.append(condata)
            except KeyError:
                # For discretization equations, which are skipped at t0
                pass
        b.cons = pyo.Reference(con_list)


def main():
    m = make_model()
    initialize_model(m)
    reshape(m)

    time = m.time
    t0 = m.time.first()
    for t in time:
        if t != t0:
            # I happen to know that conc contains the linking variables.
            # In a real implementation, we would need to find this ourselves (expensive),
            # make it a user-provided input, or make some assumption (only linking is due
            # to DAE's discretization equations).
            t_prev = time.prev(t)
            m.conc[t_prev, :].fix()
        solver.solve(m.time_block[t], tee=True)

    # Revert changes to linking variables.
    m.conc[:, :].unfix()
    m.conc[t0, :].fix()

    # Deactivate "decomposition block" as we don't want components to get
    # picked up twice when we send the full model to a solver.
    m.time_block.deactivate()

    return m


if __name__ == "__main__":
    m = main()
