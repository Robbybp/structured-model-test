import pyomo.environ as pyo
import matplotlib.pyplot as plt
import cstr_model


def plot_control(m):
    fig, ax = plt.subplots()
    time = list(m.time)
    F = list(m.flow_in[:].value)

    ax.step(time, F, linewidth=3)

    ax.tick_params(direction="in", width=3)
    for spine in ax.spines.values():
        spine.set_linewidth(3)

    plt.show()


solver = pyo.SolverFactory("ipopt")

m = cstr_model.main()

time = m.time
t0 = time.first()
t1 = time[2] # Pyomo sets are 1-indexed
comp = m.comp

m.flow_in[:].unfix()
m.flow_in[:].setlb(0.0)
m.flow_in[t0].fix()

F_tgt = 1.0
conc_setpoint = {}
conc_setpoint["A"] = F_tgt*m.conc_in[t0, "A"]/(F_tgt - m.stoich["A"])
conc_setpoint["B"] = (
        m.conc_in[t0, "B"] +
        m.stoich["B"]*conc_setpoint["A"]/F_tgt
        )

m.sp_error = pyo.Expression(time, comp, initialize={
    (t, j): m.conc[t, j] - conc_setpoint[j] for t, j in time*comp
    })

m.obj = pyo.Objective(expr=sum(m.sp_error[t, j]**2 for t in time for j in comp))

solver.solve(m, tee=True)
m.conc.pprint()
m.flow_in.pprint()

cstr_model.plot_states(m)
m.flow_in[t0].set_value(m.flow_in[t1].value)
plot_control(m)
