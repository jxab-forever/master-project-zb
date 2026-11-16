"""Additional open-source solver and validation for the paper reproduction."""
from pathlib import Path
import math
import time
import numpy as np


def check_solution(m):
    x=m.solution
    residual=0.0
    for row in m.rows:
        if not row.active: continue
        a=sum(c*x[j] for j,c in row.expression.terms.items())
        residual=max(residual,row.lower-a,a-row.upper)
    bound=max(float(np.max(np.asarray(m.lower)-x)),float(np.max(x-np.asarray(m.upper))),0.)
    integer=x[np.asarray(m.integer)]
    integ=float(np.max(abs(integer-np.round(integer)),initial=0))
    return dict(max_constraint_violation=float(max(residual,0)),max_bound_violation=bound,
                max_integrality_violation=integ)


def solve_scip(m,folder,time_limit=180,gap=.01,starts=None):
    from pyscipopt import Model,quicksum
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    scip=Model('STREAM-paper')
    scip.hideOutput()
    scip.setLogfile(str(folder/'scip.log'))
    scip.setRealParam('limits/time',float(time_limit))
    scip.setRealParam('limits/gap',float(gap))
    scip.setRealParam('numerics/feastol',1e-7)
    t=time.perf_counter()
    variables=[scip.addVar(name=n,vtype='B' if integer else 'C',
              lb=lb if math.isfinite(lb) else -scip.infinity(),
              ub=ub if math.isfinite(ub) else scip.infinity())
              for n,lb,ub,integer in zip(m.column_names,m.lower,m.upper,m.integer)]
    for row in m.rows:
        if not row.active: continue
        expression=quicksum(float(c)*variables[j] for j,c in row.expression.terms.items() if c)
        if row.lower==row.upper:
            scip.addCons(expression==row.lower,name=row.name)
        else:
            if math.isfinite(row.lower): scip.addCons(expression>=row.lower,name=row.name+'_lo')
            if math.isfinite(row.upper): scip.addCons(expression<=row.upper,name=row.name+'_hi')
    scip.setObjective(quicksum(float(c)*variables[j] for j,c in m.objective.terms.items())+m.objective.constant,'minimize')
    if starts:
        sol=scip.createSol() if len(starts)==len(variables) else scip.createPartialSol()
        for j,v in starts.items():scip.setSolVal(sol,variables[j],float(v))
        scip.addSol(sol)
    build=time.perf_counter()-t
    scip.optimize()
    if scip.getNSols()==0:
        raise RuntimeError(f'SCIP {scip.getStatus()}, no feasible solution')
    sol=scip.getBestSol()
    m.solution=np.array([scip.getSolVal(sol,v) for v in variables])
    result=dict(solver='SCIP',status=str(scip.getStatus()),objective=m.value(m.objective),
                dual_bound=float(scip.getDualbound()),mip_gap=float(scip.getGap()),
                solve_seconds=float(scip.getSolvingTime()),backend_build_seconds=build,
                variables=len(variables),constraints=sum(r.active for r in m.rows))
    result.update(check_solution(m))
    return result
