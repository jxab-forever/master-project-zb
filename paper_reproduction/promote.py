"""Adopt complete, current-model incumbents without mixing stale model bounds."""
import json
import math
import shutil
import numpy as np
from .workflow import OUT,ROOT,case_signature

PARAMETERS=['demand_growth','electricity_growth','ccs_cost_factor','target','bau']
REQUIRED=['summary.json','solution.npy','solution_variables.csv','trajectory.csv','energy.csv','states.csv','annual.csv']


def promote_case(source,dest,backup):
    new=json.loads((source/'summary.json').read_text())
    expected=case_signature(ROOT/'inputs/softwarex/f_4.xlsx',*[new[k] for k in PARAMETERS])
    if new.get('case_signature')!=expected:raise ValueError(f'Stale model/input signature: {source.name}')
    missing=[name for name in REQUIRED if not (source/name).is_file()]
    if missing:raise ValueError(f'Incomplete proposal {source.name}: {missing}')
    old=json.loads((dest/'summary.json').read_text())
    if any(new[k]!=old[k] for k in PARAMETERS):raise ValueError(f'Scenario parameters differ: {source.name} -> {dest.name}')
    if not math.isfinite(new['objective']):raise ValueError('Invalid objective')
    if max(new.get(k,float('inf')) for k in ['max_constraint_violation','max_bound_violation','max_integrality_violation'])>2e-6:
        raise ValueError(f'Unvalidated proposal: {source.name}')
    x=np.load(source/'solution.npy')
    if not np.isfinite(x).all():raise ValueError('Invalid solution vector')
    if new['objective']>old['objective']+1e-6:return False
    if not backup.exists():shutil.copytree(dest,backup)
    updated=dict(new,name=dest.name)
    if dest.name!=source.name:updated['reused_identical_case']=source.name
    # Old bounds are comparable only if their input/model identities match.
    updated['dual_bound']=max(new['dual_bound'],old['dual_bound']) if old.get('case_signature')==expected else new['dual_bound']
    updated['mip_gap']=max(0.,(updated['objective']-updated['dual_bound'])/max(abs(updated['objective']),1e-12))
    shutil.copytree(source,dest,dirs_exist_ok=True)
    (dest/'summary.json').write_text(json.dumps(updated,indent=2),encoding='utf8')
    print(f'Promoted {dest.name}: {old["objective"]:.6f} -> {new["objective"]:.6f}',flush=True)
    return True


def main():
    for summary in sorted((OUT/'refined').glob('*/summary.json')):
        name=summary.parent.name
        targets=[name]+(['scenario_4','target_09'] if name=='nominal' else [])
        for target in targets:promote_case(summary.parent,OUT/'runs'/target,OUT/'initial_milp'/target)


if __name__=='__main__':main()
