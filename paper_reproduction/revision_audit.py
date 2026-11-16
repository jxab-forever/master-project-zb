"""Check solution vectors directly against independently exported Julia MPS files.

This checker does not construct equations with the Python model builder and does
not trust residuals recorded in old summary files. Only the cumulative carbon
row is adjusted for target sweeps or disabled for BAU, exactly as in Julia.
"""
from pathlib import Path
import csv
import hashlib
import json
import os
import numpy as np
import pandas as pd
import highspy
from scipy.sparse import csc_matrix
from validation.compare_julia import canonical
from .workflow import ROOT,OUT,FACTORS,TARGETS

_MODELS={}


def reference_model(path):
    path=Path(path).resolve()
    if path in _MODELS:return _MODELS[path]
    h=highspy.Highs();h.setOptionValue('output_flag',False)
    old=os.getcwd()
    try:
        os.chdir(ROOT)
        status=h.readModel(str(path.relative_to(ROOT)))
    finally:os.chdir(old)
    if status==highspy.HighsStatus.kError:raise RuntimeError(f'Cannot read Julia MPS: {path}')
    lp=h.getLp()
    matrix=csc_matrix((lp.a_matrix_.value_,lp.a_matrix_.index_,lp.a_matrix_.start_),shape=(lp.num_row_,lp.num_col_)).tocsr()
    scalar={}
    parameter_file=path.with_name(path.name.replace('_julia.mps','_julia_parameters.csv'))
    for row in csv.DictReader(parameter_file.open(encoding='utf-8-sig',newline='')):
        if row['index']=='0':scalar[row['parameter']]=float(row['value'])
    carbon=[i for i,n in enumerate(lp.row_names_) if canonical(n)=='ag_co2_l_link_i_']
    if len(carbon)!=1:raise RuntimeError('Expected exactly one cumulative carbon row')
    model=dict(lp=lp,matrix=matrix,carbon=carbon[0],scalar=scalar,
        names=[canonical(n) for n in lp.col_names_],sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    _MODELS[path]=model
    return model


def validate_case(folder,reference,target,bau,expected_parameters):
    folder=Path(folder);model=reference_model(reference);lp=model['lp'];scalar=model['scalar']
    summary=json.loads((folder/'summary.json').read_text())
    df=pd.read_csv(folder/'solution_variables.csv');stored=np.load(folder/'solution.npy')
    names=[canonical(n) for n in df.variable]
    if len(names)!=len(set(names)) or set(names)!=set(model['names']):
        return dict(case=folder.name,passed=False,error='Variable names do not match Julia model')
    lookup={n:i for i,n in enumerate(names)}
    values=df.value.to_numpy(float)
    if stored.shape!=values.shape or not np.isfinite(values).all():
        return dict(case=folder.name,passed=False,error='Invalid or incomplete solution vector')
    x=values[[lookup[n] for n in model['names']]]
    lower=np.array(lp.row_lower_);upper=np.array(lp.row_upper_)
    c=model['carbon'];years=scalar['n_periods']*scalar['n_subperiods']*scalar['yr_subperiod']
    upper[c]=scalar['initial_co2']*years*(1-float(target)/2)
    if bau:lower[c]=-np.inf;upper[c]=np.inf
    activity=model['matrix']@x
    row_violation=np.maximum(np.maximum(lower-activity,activity-upper),0.)
    maxrow=float(np.max(row_violation,initial=0))
    maxbound=max(float(np.max(np.array(lp.col_lower_)-x,initial=0)),float(np.max(x-np.array(lp.col_upper_),initial=0)))
    integer=np.array([int(v)!=0 for v in lp.integrality_],bool)
    maxinteger=float(np.max(abs(x[integer]-np.round(x[integer])),initial=0))
    computed_objective=float(np.dot(lp.col_cost_,x)+lp.offset_)
    cost_difference=abs(computed_objective-summary['objective'])
    vector_difference=float(np.max(abs(stored-values),initial=0))
    params=(summary['demand_growth'],summary['electricity_growth'],summary['ccs_cost_factor'])
    metadata_ok=bool(np.allclose(params,expected_parameters,rtol=0,atol=1e-14) and abs(summary['target']-target)<1e-14 and bool(summary['bau'])==bool(bau))
    computed_co2=float(activity[c]/scalar['sf_em'])
    co2_difference=abs(computed_co2-summary['total_co2'])
    gap=max(0.,(computed_objective-summary['dual_bound'])/max(abs(computed_objective),1e-12))
    checks=dict(original_julia_rows=maxrow<2e-6,variable_bounds=maxbound<2e-6,
        binary_integrality=maxinteger<1e-7,objective_accounting=cost_difference<2e-5,
        saved_vector_consistency=vector_difference<1e-10,scenario_metadata=metadata_ok,
        cumulative_emissions=co2_difference<2e-5,lower_bound_not_above_incumbent=summary['dual_bound']<=computed_objective+2e-5,
        gap_accounting=abs(gap-summary['mip_gap'])<1e-6)
    checks={k:bool(v) for k,v in checks.items()}
    worst=int(np.argmax(row_violation))
    return dict(case=folder.name,passed=all(checks.values()),checks=checks,
        max_constraint_violation=maxrow,max_bound_violation=maxbound,max_integrality_violation=maxinteger,
        objective=computed_objective,objective_difference=cost_difference,total_co2=computed_co2,
        emissions_difference=co2_difference,saved_vector_difference=vector_difference,
        mip_gap=gap,worst_row=lp.row_names_[worst],reference=Path(reference).name,reference_sha256=model['sha256'])


def audit_all():
    references=ROOT/'validation/revision_20261110'
    cases=[('nominal',4,.5,False),('bau_nominal',4,.5,True)]
    cases += [(f'scenario_{n}{suffix}',n,.5,bau) for n in range(1,9) for suffix,bau in [('',False),('_bau',True)]]
    cases += [(f'target_{i:02d}',4,float(t),False) for i,t in enumerate(TARGETS)]
    results=[]
    for name,n,target,bau in cases:
        result=validate_case(OUT/'runs'/name,references/f'scenario_{n}_julia.mps',target,bau,FACTORS[n-1])
        results.append(result);print(name,'PASS' if result['passed'] else 'FAIL',flush=True)
    report=dict(passed=all(r['passed'] for r in results),cases=len(results),
        method='Independent Julia MPS matrices; Python equation builder not used by checker',
        max_constraint_violation=max(r.get('max_constraint_violation',float('inf')) for r in results),results=results)
    (OUT/'independent_validation.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    if not report['passed']:raise RuntimeError('Independent Julia solution audit failed')
    return report


if __name__=='__main__':audit_all()
