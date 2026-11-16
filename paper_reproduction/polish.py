"""Second-pass joint timing search using immutable first-pass neighboring plans."""
import argparse,json
import numpy as np
from .workflow import ROOT,OUT,scenario_params,sets,build_softwarex,export_data
from .joint_search import search_joint_plans
from .solver import check_solution


def polish(name,seconds=90):
    folder=OUT/'refined'/name
    summary=json.loads((folder/'summary.json').read_text())
    p=scenario_params(ROOT/'inputs/softwarex/f_4.xlsx',summary['demand_growth'],summary['electricity_growth'],summary['ccs_cost_factor'],summary['target'])
    s=sets(p);m=build_softwarex(p,s)
    if summary['bau']:m.remove_group('ag_co2_l_link_i_')
    x=np.load(folder/'solution.npy');initial=summary['objective']
    bank=OUT/'revision_20261110/first_pass'
    if name.startswith('target_'):
        i=int(name.split('_')[1]);neighbors=[f'target_{j:02d}' for j in [i+1,i+2,i-1] if 0<=j<=8]
    else:
        i=int(name.split('_')[1]);neighbors=[f'scenario_{i+1 if i%2 else i-1}']
        if i<=3:neighbors+=['nominal','scenario_3']
    seeds=[np.load(bank/n/'solution.npy') for n in neighbors if n!=name and (bank/n/'solution.npy').exists()]
    best,stats=search_joint_plans(m,p,s,x,seconds=seconds,strategies=('timing','coordinate'),seed_solutions=seeds)
    value=m.value(m.objective)
    if value>initial+1e-6:raise RuntimeError('Polishing worsened the validated incumbent')
    passes=summary.get('refinement_passes',[summary['refinement']])
    summary.update(check_solution(m))
    summary.update(objective=value,primal_bound=value,mip_gap=max(0.,(value-summary['dual_bound'])/abs(value)),
        refinement=stats,refinement_passes=passes+[stats],search_version='joint_timing_v2_20261110')
    export_data(m,p,s,folder,summary)
    print(f'POLISHED {name}: {initial:.6f} -> {value:.6f}',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--group',choices=['targets','factorial'],required=True);p.add_argument('--seconds',type=float,default=90);args=p.parse_args()
    names=[f'target_{i:02d}' for i in range(9)] if args.group=='targets' else [f'scenario_{i}' for i in [1,2,3,5,6,7,8]]
    for name in names:polish(name,args.seconds)
