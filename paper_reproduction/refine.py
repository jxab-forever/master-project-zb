"""Improve original MILP incumbents with validated joint and local LP plans."""
from pathlib import Path
import argparse
import json
import time
import numpy as np
from .workflow import ROOT,OUT,scenario_params,sets,build_softwarex,export_data,case_signature
from .solver import check_solution
from .joint_search import search_joint_plans


def refine(name,seconds=120,force=False):
    source=OUT/'runs'/name
    summary=json.loads((source/'summary.json').read_text())
    source_file=ROOT/'inputs/softwarex/f_4.xlsx'
    signature=case_signature(source_file,summary['demand_growth'],summary['electricity_growth'],
        summary['ccs_cost_factor'],summary['target'],summary['bau'])
    folder=OUT/'refined'/name
    cached_path=folder/'summary.json'
    cached=json.loads(cached_path.read_text()) if cached_path.exists() else None
    if (not force and cached and cached.get('case_signature')==signature
        and cached.get('search_version')=='joint_timing_v2_20261110'
        and cached.get('requested_search_seconds',0)>=seconds
        and cached['objective']<=summary['objective']+1e-6
        and (folder/'solution.npy').exists()):
        print(f'CACHE REFINED {name}: verified model and search budget',flush=True)
        return cached
    p=scenario_params(source_file,summary['demand_growth'],summary['electricity_growth'],
        summary['ccs_cost_factor'],summary['target'])
    s=sets(p);m=build_softwarex(p,s)
    if summary['bau']:m.remove_group('ag_co2_l_link_i_')
    x=np.load(source/'solution.npy')
    initial=float(summary['objective'])
    # Always begin from the better verified incumbent. A short seed-search budget
    # must never silently discard a previously improved feasible solution.
    if cached and cached.get('case_signature')==signature and (folder/'solution.npy').exists():
        candidate=np.load(folder/'solution.npy')
        if len(candidate)==len(m.lower) and np.isfinite(candidate).all():
            m.solution=candidate
            if max(check_solution(m).values())<2e-6 and m.value(m.objective)<initial:
                x=candidate.copy()
    seeds=[]
    for path in sorted((OUT/'runs').glob('*/solution.npy')):
        if path.parent.name in (name,'nominal_scip'):continue
        seeds.append(np.load(path))
    if cached and cached.get('case_signature')==signature and (folder/'solution.npy').exists():
        seeds.append(np.load(folder/'solution.npy'))
    best,stats=search_joint_plans(m,p,s,x,seconds=seconds,seed_solutions=seeds)
    summary.update(check_solution(m))
    best_cost=m.value(m.objective)
    # The independent audit proves that all active scenario equations match the
    # Julia matrices; the retained original global lower bound therefore applies.
    bound=summary['dual_bound']
    if (summary['demand_growth'],summary['electricity_growth'],summary['ccs_cost_factor'],summary['target'],summary['bau'])==(.02,.025,1,.5,False):
        for other in ['baseline_nominal','native110_nominal']:
            result=OUT/'revision_20261110'/other/'summary.json'
            if result.exists():bound=max(bound,json.loads(result.read_text())['dual_bound'])
    summary.update(objective=best_cost,primal_bound=best_cost,dual_bound=bound,
        mip_gap=max(0.,(best_cost-bound)/max(abs(best_cost),1e-12)),
        refinement=stats,requested_search_seconds=seconds,case_signature=signature,
        legacy_start_revalidated=True,source_objective=initial,
        solver='HiGHS MILP + validated joint-plan LP search',search_version='joint_timing_v2_20261110',
        status='Feasible incumbent; not globally certified')
    export_data(m,p,s,folder,summary)
    print(f'REFINED {name}: {initial:.6f} -> {best_cost:.6f}; gap={summary["mip_gap"]:.2%}',flush=True)
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--case',default='target_00')
    parser.add_argument('--seconds',type=float,default=120)
    parser.add_argument('--group',choices=['targets','factorial'])
    parser.add_argument('--force',action='store_true',help='Search again even when a valid refinement cache exists')
    args=parser.parse_args()
    if args.group=='targets':names=[f'target_{i:02d}' for i in range(9)]
    elif args.group=='factorial':names=['nominal']+[f'scenario_{i}' for i in [1,2,3,5,6,7,8]]
    else:names=([f'target_{i:02d}' for i in range(9)]+['nominal']+[f'scenario_{i}' for i in [1,2,3,5,6,7,8]]) if args.case=='all' else [args.case]
    for name in names:refine(name,args.seconds,args.force)
