import pyomo.environ as pyo
import pyomo.dae as dae


solver = pyo.SolverFactory("ipopt")


def _block_rule(b, t):
    conc = b.parent_block().conc
    dcdt = b.parent_block().dcdt
    comp = b.parent_block().comp
    stoich = b.parent_block().stoich
    k_rxn = b.parent_block().k_rxn



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

    m.flow_eqn = pyo.Constraint(expr=m.flow_in[t] == m.flow_out[t])

    m.conc_in = pyo.Var(time, comp)
    m.conc_out = pyo.Var(time, comp)
    m.conc_out_eqn = pyo.Constraint(time, comp, rule={
        j: conc[t, j] - m.conc_out[t, j] == 0 for j in m.comp for t in time
        })

    m.rate_gen = pyo.Var(time, comp)
    m.rate_eqn = pyo.Constraint(time, comp, rule={
        j: m.rate_gen[t, j] - m.stoich[j]*m.k_rxn*m.conc[t, "A"] == 0
        for j in m.comp for t in time
        })

    m.conc_diff_eqn = pyo.Constraint(time, comp, rule={
        j: dcdt[t, j] -
           (
               m.flow_in[t]*m.conc_in[t, j] -
               m.flow_out[t]*m.conc_out[t, j] +
               m.rate_gen[t, j]
               ) == 0
        for j in m.comp
        })


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
