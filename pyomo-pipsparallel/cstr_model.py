import pyomo.environ as pyo
import pyomo.dae as dae


solver = pyo.SolverFactory("ipopt")


def _block_rule(b, t):
    conc = b.parent_block().conc
    dcdt = b.parent_block().dcdt
    comp = b.parent_block().comp
    stoich = b.parent_block().stoich
    k_rxn = b.parent_block().k_rxn

    b.flow_in = pyo.Var()
    b.flow_out = pyo.Var()
    b.flow_eqn = pyo.Constraint(expr=b.flow_in == b.flow_out)

    b.conc_in = pyo.Var(comp)
    b.conc_out = pyo.Var(comp)
    b.conc_out_eqn = pyo.Constraint(comp, rule={
        j: conc[t, j] - b.conc_out[j] == 0 for j in comp
        })

    b.rate_gen = pyo.Var(comp)
    b.rate_eqn = pyo.Constraint(comp, rule={
        j: b.rate_gen[j] - stoich[j]*k_rxn*conc[t, "A"] == 0 for j in comp
        })

    b.conc_diff_eqn = pyo.Constraint(comp, rule={
        j: dcdt[t, j] -
           (
               b.flow_in*b.conc_in[j] -
               b.flow_out*b.conc_out[j] +
               b.rate_gen[j]
               ) == 0
        for j in comp
        })


def make_model(horizon=10.0):
    m = pyo.ConcreteModel()
    m.comp = pyo.Set(initialize=["A", "B"])
    m.time = dae.ContinuousSet(initialize=[0, horizon])

    m.stoich = pyo.Param(m.comp, initialize={"A": -1, "B": 1}, mutable=True)
    m.k_rxn = pyo.Param(initialize=1.0, mutable=True)

    m.conc = pyo.Var(m.time, m.comp)
    m.dcdt = dae.DerivativeVar(m.conc, wrt=m.time)

    m.time_block = pyo.Block(m.time, rule=_block_rule)

    return m


def initialize_model(m, ntfe=10):
    disc = pyo.TransformationFactory("dae.finite_difference")
    disc.apply_to(m, wrt=m.time, nfe=ntfe, scheme="BACKWARD")

    t0 = m.time.first()

    m.time_block[:].conc_in["A"].fix(5.0)
    m.time_block[:].conc_in["B"].fix(0.01)

    m.time_block[:].flow_in.fix(1.0)
    m.time_block[t0].flow_in.fix(0.1)

    m.conc[t0, "A"].fix(1.0)
    m.conc[t0, "B"].fix(0.0)


def main():
    m = make_model()
    initialize_model(m)
    solver.solve(m, tee=True)
    return m


if __name__ == "__main__":
    m = main()
    import pdb; pdb.set_trace()
