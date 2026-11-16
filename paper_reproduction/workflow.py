"""Rebuild and solve the paper's BAU, target sweep, and factorial scenarios.

Copyright (C) 2023 UChicago Argonne, LLC, for the underlying STREAM equations.
Scenario transformations follow gen_data.jl in the paper's pinned repository.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import json
import math
import time
import shutil
import hashlib
import inspect
import numpy as np
import pandas as pd
from stream_python.softwarex_data import read_params,sets,JArray
from stream_python.softwarex_model import build_softwarex,gen_lp_window
from .solver import check_solution

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'paper_reproduction'
FACTORS=[(.02,.05,2),(.02,.05,1),(.02,.025,2),(.02,.025,1),
         (0,.05,2),(0,.05,1),(0,.025,2),(0,.025,1)]
TARGETS=np.linspace(0,.5,10)
REFRESH=False


def scenario_params(source,demand_growth,electricity_growth,ccs_cost,target):
    p=read_params(source)
    if (demand_growth,electricity_growth,ccs_cost,target)==(.02,.025,1,.5):return p
    put=lambda name,a:setattr(p,name,JArray(np.asarray(a,dtype=float)))
    years=np.arange(p.n_periods*p.n_subperiods)*p.yr_subperiod
    initial=float(p.demand[1,1])
    final=initial*(1+demand_growth)**int(years[-1])
    put('demand',np.linspace(initial,final,len(years)).reshape(p.n_periods,p.n_subperiods))
    put('elec_budget',p.elec_budget[1,1]*(1+electricity_growth)**years.reshape(p.n_periods,p.n_subperiods))
    put('co2_budget',np.full((p.n_periods,p.n_subperiods),p.initial_co2*(1-target/2)*len(years)*p.yr_subperiod))
    # Recompute physical capacity bounds, including those generated from demand.
    expansion=math.ceil((final-initial)/initial)
    new_max=math.ceil(final/initial)
    old_cpb=p.r_cpb_ub.data.copy();old_nc0=p.n_c0_bM.data.copy()
    put('e_C',np.floor(p.c0.data/p.sf_cap)*p.sf_cap*expansion/p.x_ub)
    put('e_c_ub',p.e_C.data*p.x_ub*1.1)
    put('r_cpb_ub',(p.e_c_ub.data+p.c0.data)*1.1)
    rratio=p.r_cpb_ub.data/old_cpb
    put('n_c0_bM',p.c0.data*new_max)
    nratio=p.n_c0_bM.data/old_nc0
    for prefix,ratio,names in [('r',rratio,['cp_ub','eh_ub','u_ub','cpe_ub','fu_e_ub','u_fu_e_ub',
                 'ep0_ub','ep1ge_ub','ep1gce_ub','ep1gcs_ub','ups_e_mt_in_ub','x_in_ub','x_out_ub','cvonm_ub']),
                ('n',nratio,['cp_bM','eh_ub','u_ub','cpe_ub','fu_e_ub','u_fu_e_ub','ep0_bM',
                 'ep1ge_bM','ep1gce_bM','ep1gcs_bM','ups_e_mt_in_ub','x_in_ub','x_out_ub','cvonm_bM','fstck_ub'])]:
        for suffix in names:
            name=prefix+'_'+suffix;a=getattr(p,name).data
            put(name,a*ratio.reshape((-1,)+(1,)*(a.ndim-1)))
        for suffix in ['ehf_ub','u_ehf_ub']:
            name=prefix+'_'+suffix;a=getattr(p,name).data
            put(name,a*ratio.reshape((1,-1)+(1,)*(a.ndim-2)))
    put('r_e_c_ub',p.e_C.data*p.x_ub)
    # Capital multipliers also change fixed O&M and all dependent cash bounds.
    p.r_loanFact.data[:,:,:,4]*=ccs_cost
    p.n_loanFact.data[:,4]*=ccs_cost
    put('r_c_fOnm',p.r_loanFact.data[0,0]*.02)
    put('n_c_fOnm',p.n_loanFact.data*.02)
    ann=float(p.e_Ann[1])
    put('e_l_ub',p.e_c_ub.data*p.e_loanFact.data)
    for name in ['e_ladd_ub','e_loan_ub']:put(name,p.e_l_ub.data)
    for name in ['e_ann_ub','e_pay_ub']:put(name,p.e_l_ub.data*ann)
    put('r_l0_ub',p.c0.data[:,None]*p.r_loanFact.data[0,0]*1.1)
    put('r_le_ub',(p.e_C.data*p.x_ub)[:,None]*p.r_loanFact.data[0,0]*1.1)
    rloan=p.r_cpb_ub.data*np.max(p.r_loanFact.data[0,0],axis=1)
    for name in ['r_l0add_bM','r_leadd_bM']:put(name,rloan)
    put('r_loan_ub',rloan+p.r_loan0.data)
    for name in ['r_ann0_bM','r_pay0_bM']:put(name,np.max(p.r_l0_ub.data,axis=1)*ann)
    for name in ['r_anne_bM','r_paye_bM']:put(name,np.max(p.r_le_ub.data,axis=1)*ann)
    put('r_cfonm_ub',p.r_cpb_ub.data*np.max(p.r_c_fOnm.data,axis=1))
    for name in ['t_loan_bM','t_ret_c_bM']:put(name,p.r_loan_ub.data+p.e_loan_ub.data)
    put('o_pay_bM',rloan*ann+p.r_anne_bM.data+p.e_ann_ub.data)
    put('o_cfonm_bM',p.r_cfonm_ub.data)
    put('o_cvonm_bM',p.r_cvonm_ub.data)
    put('n_l_bM',p.n_c0_bM.data*np.max(p.n_loanFact.data,axis=1))
    for name in ['n_ladd_bM','n_loan_bM']:put(name,p.n_l_bM.data)
    for name in ['n_ann_bM','n_pay_bM']:put(name,p.n_l_bM.data*ann)
    put('n_cfonm_bM',p.n_c0_bM.data*np.max(p.n_c_fOnm.data))
    for target_name,source_name in [('o_cp_ub','r_cp_ub'),('o_u_ub','r_u_ub'),('o_ep0_bM','r_ep0_ub'),
        ('o_ep1ge_bM','r_ep1ge_ub'),('o_ep1gce_bM','r_ep1gce_ub'),('o_ep1gcs_bM','r_ep1gcs_ub'),
        ('o_ups_e_mt_in_ub','r_ups_e_mt_in_ub'),('o_x_in_ub','r_x_in_ub'),('o_x_out_ub','r_x_out_ub')]:
        put(target_name,getattr(p,source_name).data)
    put('o_ehf_ub',p.r_ehf_ub.data+p.r_u_ehf_ub.data)
    # gen_data.jl:1224, including this currently inactive feedstock bound.
    put('o_fstck_ub',np.max(p.r_fstck_ub.data,axis=1)+np.max(p.n_fstck_ub.data,axis=1))
    return p


def plan_values(m,retire=None,retro=None,new=None,expand=None):
    retire=retire or {};retro=retro or {};new=new or {};expand=expand or {}
    starts={}
    for family in ['y_o','y_e','y_r','y_n']:
        for key,var in m.variables[family].items():
            i,j,l=key[:3];t=(i-1)*4+j
            if family=='y_o':v=t<retire.get(l,99)
            elif family=='y_e':v=t>=expand.get(l,99)
            elif family=='y_r':
                tr,k=retro.get(l,(99,1));v=key[3]==(k if t>=tr else 1)
            else:
                tr,k=new.get(l,(99,1));v=key[3]==(k if t>=tr else 1)
            starts[var.index]=float(v)
    return starts


def source_binary_values(m):
    starts={}
    for family,f in [('y_o','yostate.csv'),('y_e','yestate.csv'),('y_r','yrstate.csv'),('y_n','ynstate.csv')]:
        df=pd.read_csv(ROOT/'inputs/softwarex/warmstart'/f)
        for row in df.itertuples(index=False,name=None):
            key=tuple(int(x) for x in row[:-1])
            if key in m.variables[family]:starts[m.variables[family][key].index]=float(row[-1])
    return starts


def fixed_lp(m,starts):
    lower,upper,integer=m.lower,m.upper,m.integer
    m.lower=lower.copy();m.upper=upper.copy();m.integer=[False]*len(integer)
    for k,v in starts.items():m.lower[k]=m.upper[k]=round(v)
    try:
        r=m.solve(time_limit=15,mip_gap=0,threads=1)
        return m.solution.copy(),r['objective']
    except RuntimeError:return None,float('inf')
    finally:m.lower,m.upper,m.integer=lower,upper,integer


def export_data(m,p,s,folder,summary):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    value=lambda family,key: max(0.,m.value(m.variables[family][key])) if key in m.variables[family] else 0.
    rows=[];fuel=[];events=[];annual=[]
    for i in s.P:
        for j in s.P2:
            t=(i-1)*p.n_subperiods+j
            for l in s.L:
                online=round(value('y_o',(i,j,l)))
                for form,ks,fltr in [('existing',s.Kr,p.r_filter),('new',s.Kn,p.n_filter)]:
                    for k in ks:
                        if not fltr[l,k]:continue
                        if form=='existing':
                            cap=value('r_cpb_d_',(i,j,l,k))*online/p.sf_cap
                            production=value('r_cp_d_',(i,j,l,k,p.key_node))*online/p.sf_cap
                            emit=value('r_ep1ge_d_',(i,j,l,k))*online/p.sf_em
                            captured=(value('r_ep1gce_d_',(i,j,l,k))+value('r_ep1gcs_d_',(i,j,l,k)))*online/p.sf_em
                            energy=sum(value('r_u_d_',(i,j,l,k,n)) for n in s.Nd if p.nd_en_fltr[n])*online/p.sf_elec
                        else:
                            cap=value('n_c0_d_',(i,j,l,k))/p.sf_cap
                            production=value('n_cp_d_',(i,j,l,k,p.key_node))/p.sf_cap
                            emit=value('n_ep1ge_d_',(i,j,l,k))/p.sf_em
                            captured=(value('n_ep1gce_d_',(i,j,l,k))+value('n_ep1gcs_d_',(i,j,l,k)))/p.sf_em
                            energy=sum(value('n_u_d_',(i,j,l,k,n)) for n in s.Nd if p.nd_en_fltr[n])/p.sf_elec
                        rows.append(dict(subperiod=t,location=l,form=form,technology=k,
                            installed_capacity=cap,production=production,direct_emissions=emit,captured=captured))
                        fuel.append(dict(subperiod=t,location=l,form=form,technology=k,fuel='Electricity',energy=energy,
                            paper_energy=energy,process_energy=energy,onsite_generation_energy=0.))
                        for f,label in [(1,'Coal'),(2,'Fuel oil'),(3,'Natural gas')]:
                            prefix='r' if form=='existing' else 'n'
                            scale=(online if form=='existing' else 1)/p.sf_heat
                            process=sum(value(prefix+'_ehf_d_',(i,j,l,k,f,n))
                                for n in s.Nd if p.nd_en_fltr[n])*scale
                            onsite=sum(value(prefix+'_u_ehf_d_',(i,j,l,k,f,n))
                                for n in s.Nd if p.nd_en_fltr[n])*scale
                            # Published write_fuel_res_v1 excludes onsite-generation
                            # fuel. Retain the full scope separately for accounting.
                            fuel.append(dict(subperiod=t,location=l,form=form,technology=k,fuel=label,
                                energy=process+onsite,paper_energy=process,
                                process_energy=process,onsite_generation_energy=onsite))
                ex_tech=max(s.Kr,key=lambda k:value('y_r',(i,j,l,k)))
                new_tech=max(s.Kn,key=lambda k:value('y_n',(i,j,l,k)))
                events.append(dict(subperiod=t,location=l,online=online,
                    expansion=value('e_c',(i,j,l))*online/p.sf_cap,
                    existing_technology=ex_tech,new_technology=new_tech))
            a=dict(subperiod=t,demand=p.demand[i,j]/p.sf_cap,
                electricity_limit=p.elec_budget[i,j]/p.sf_elec,capture_limit=p.initial_co2*.9/p.sf_em)
            a['budget_emissions']=sum(value(f,(i,j,l)) for l in s.L for f in
                ['o_ep1ge','n_ep1ge','o_ups_e_mt_in','n_ups_e_mt_in'])/p.sf_em
            annual.append(a)
    d=pd.DataFrame(rows);e=pd.DataFrame(fuel);a=pd.DataFrame(annual)
    summary.update(total_co2=float(a.budget_emissions.sum()*p.yr_subperiod),
        min_demand_margin=float((d.groupby('subperiod').production.sum().to_numpy()-a.demand).min()))
    loan,_=gen_lp_window(p,s);i=p.n_periods;j=p.n_subperiods
    terminal=.001*p.discount[i,j]*(sum(value('o_loan_last',(i,l,1)) for l in s.L)+
        sum(value('n_loan',(i,j,i0,j0,l)) for l in s.L for i0 in s.P for j0 in s.P2
            if (i,j) in loan[i0,j0] and (i,j)!=(i0,j0)))
    summary.update(terminal_penalty=float(terminal),npv=float(summary['objective']-terminal))
    for label,df in [('trajectory',d),('energy',e),('states',pd.DataFrame(events)),('annual',a)]:df.to_csv(folder/f'{label}.csv',index=False)
    np.save(folder/'solution.npy',m.solution)
    pd.DataFrame({'variable':m.column_names,'value':m.solution}).to_csv(folder/'solution_variables.csv',index=False)
    (folder/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf8')
    return summary


def case_signature(source,demand,electricity,ccs,target,bau=False):
    """Identity of the mathematical inputs; report formatting is not an input."""
    files=[Path(source),ROOT/'stream_python/softwarex_data.py',
        ROOT/'stream_python/softwarex_model.py',ROOT/'stream_python/linear.py']
    payload=dict(schema=3,parameters=[float(demand),float(electricity),float(ccs),float(target),bool(bau)],
        scenario_constructor_sha256=hashlib.sha256(inspect.getsource(scenario_params).encode()).hexdigest(),
        sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    return hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()


def completed_case(folder,signature,seconds=None,gap=.01):
    required=['summary.json','solution.npy','solution_variables.csv','trajectory.csv','energy.csv','states.csv','annual.csv']
    if not all((folder/name).is_file() for name in required):return None
    try:summary=json.loads((folder/'summary.json').read_text())
    except (OSError,ValueError):return None
    if summary.get('case_signature')!=signature:return None
    if seconds is not None and summary.get('mip_gap',1)>gap and summary.get('time_limit',0)<seconds:return None
    return summary


def solve_case(name,demand,electricity,ccs,target,bau=False,seconds=300,mip_gap=.01,threads=2):
    folder=OUT/'runs'/name
    source=ROOT/'inputs/softwarex/f_4.xlsx'
    signature=case_signature(source,demand,electricity,ccs,target,bau)
    cached=completed_case(folder,signature,seconds,mip_gap)
    if cached is not None and not REFRESH:
        print(f'CACHE {name}: verified input/model signature',flush=True)
        return cached
    if name not in ['nominal','bau_nominal'] and (demand,electricity,ccs,target)==(.02,.025,1,.5):
        equivalent=OUT/'runs'/('bau_nominal' if bau else 'nominal')
        if completed_case(equivalent,signature,seconds,mip_gap) is not None:
            shutil.copytree(equivalent,folder,dirs_exist_ok=True)
            summary=json.loads((folder/'summary.json').read_text());summary['name']=name
            summary['reused_identical_case']=equivalent.name
            (folder/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf8')
            return summary
    print(f'RUN {name}: D={demand} E={electricity} C={ccs} target={target} BAU={bau}',flush=True)
    p=scenario_params(ROOT/'inputs/softwarex/f_4.xlsx',demand,electricity,ccs,target)
    s=sets(p);m=build_softwarex(p,s)
    if bau:m.remove_group('ag_co2_l_link_i_')
    candidates=[source_binary_values(m)]
    expand={l:1 for l in s.L} if demand else {}
    candidates.append(plan_values(m,expand=expand))
    # A visible plan in the published figure is an initial candidate only;
    # all quantities below are re-optimized, and the final MILP remains free.
    if demand:
        candidates.append(plan_values(m,retire={1:2,2:5,5:8},retro={3:(6,4),4:(5,5)},
            new={1:(4,5),2:(5,5),5:(8,4)},expand={2:2,3:7,4:6}))
    for onset in ([3,5,7] if demand else [3,5,7]):
        candidates.append(plan_values(m,retro={1:(onset,5)},expand=expand))
    if demand:
        for onset in [3,5]:
            candidates.append(plan_values(m,retro={1:(onset,5),4:(onset,5),3:(6,4)},expand=expand))
    else:
        # Without demand growth the generator permits no incumbent expansion.
        # Replacement (new CCS) can remain feasible where retrofit-only seeds fail.
        for onset in [3,4,5,6]:
            candidates.append(plan_values(m,retire={1:onset},new={1:(onset,5)}))
        candidates.append(plan_values(m,retire={1:6,2:6},new={1:(6,5),2:(6,5)}))
    # Neighboring budget/scenario runs offer valid starting patterns after LP repair.
    for other in sorted((OUT/'runs').glob('*/solution.npy')):
        if other.parent.name=='nominal_scip':continue
        x=np.load(other)
        if len(x)==len(m.integer):candidates.append({i:round(x[i]) for i,b in enumerate(m.integer) if b})
    seen=set();best=None;cost=float('inf');tested=0
    for candidate in candidates:
        key=tuple(candidate.values())
        if key in seen:continue
        seen.add(key);tested+=1
        x,c=fixed_lp(m,candidate)
        if x is not None and c<cost:best,cost=x,c
    if best is None:raise RuntimeError(f'No feasible initial plan for {name}')
    print(f'{name}: {tested} LP-repaired candidates; incumbent={cost:.6f}',flush=True)
    folder.mkdir(parents=True,exist_ok=True)
    # Native solver file paths are ASCII relative names, for Windows portability.
    summary=m.solve(time_limit=seconds,mip_gap=mip_gap,threads=threads,
        log_file=folder/'highs.log',starts={i:float(v) for i,v in enumerate(best)})
    if summary['objective']>cost+1e-5:
        raise RuntimeError('Solver lost its validated initial incumbent')
    summary.update(name=name,demand_growth=demand,electricity_growth=electricity,
        ccs_cost_factor=ccs,target=target,bau=bau,initial_plans_checked=tested,
        budget=None if bau else p.co2_budget[1,1]/p.sf_em,
        solver='HiGHS',time_limit=seconds,case_signature=signature,requested_gap=mip_gap,threads=threads)
    print(f'DONE {name}: {summary["objective"]:.6f}; gap={summary["mip_gap"]:.2%}',flush=True)
    return export_data(m,p,s,folder,summary)


def main():
    global REFRESH, OUT
    parser=argparse.ArgumentParser();parser.add_argument('--seconds',type=float,default=300)
    parser.add_argument('--only',choices=['all','nominal','factorial','pareto'],default='all')
    parser.add_argument('--refresh',action='store_true',help='Re-solve completed checkpoints')
    parser.add_argument('--output',type=Path,default=OUT,help='Directory for runs and summaries')
    parser.add_argument('--gap',type=float,default=.01)
    parser.add_argument('--threads',type=int,default=2)
    args=parser.parse_args()
    REFRESH=args.refresh
    OUT=args.output.resolve()
    if args.only in ['all','nominal']:
        solve_case('nominal',.02,.025,1,.5,seconds=args.seconds,mip_gap=args.gap,threads=args.threads)
        solve_case('bau_nominal',.02,.025,1,.5,bau=True,seconds=args.seconds,mip_gap=args.gap,threads=args.threads)
    if args.only in ['all','factorial']:
        for n,(d,e,c) in enumerate(FACTORS,1):
            solve_case(f'scenario_{n}',d,e,c,.5,seconds=args.seconds,mip_gap=args.gap,threads=args.threads)
            solve_case(f'scenario_{n}_bau',d,e,c,.5,bau=True,seconds=args.seconds,mip_gap=args.gap,threads=args.threads)
    if args.only in ['all','pareto']:
        for n,target in enumerate(TARGETS):
            solve_case(f'target_{n:02d}',.02,.025,1,float(target),seconds=args.seconds,mip_gap=args.gap,threads=args.threads)
    summaries=[json.loads(p.read_text()) for p in (OUT/'runs').glob('*/summary.json')]
    pd.DataFrame(summaries).to_csv(OUT/'scenario_summary.csv',index=False)


if __name__=='__main__':main()
