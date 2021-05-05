import pyomo.environ as pyo
import cstr_model


solver = pyo.SolverFactory("ipopt")

m = cstr_model.main()

time = m.time
t0 = time.first()
comp = m.comp

m.time_block[:].flow_in.unfix()
m.time_block[:].flow_in.setlb(0.0)
m.time_block[t0].flow_in.fix()

F_tgt = 1.0
conc_setpoint = {}
conc_setpoint["A"] = F_tgt*m.time_block[t0].conc_in["A"]/(F_tgt - m.stoich["A"])
conc_setpoint["B"] = (
        m.time_block[t0].conc_in["B"] +
        m.stoich["B"]*conc_setpoint["A"]/F_tgt
        )

m.sp_error = pyo.Expression(time, comp, initialize={
    (t, j): m.conc[t, j] - conc_setpoint[j] for t, j in time*comp
    })

m.obj = pyo.Objective(expr=sum(m.sp_error[t, j]**2 for t in time for j in comp))

solver.solve(m, tee=True)

m.conc.pprint()
m.time_block[:].flow_in.pprint()

import pdb; pdb.set_trace()
