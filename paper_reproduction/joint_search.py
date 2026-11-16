"""Validated joint-plan primal search for the unchanged SoftwareX MILP.

Every candidate fixes ALL binary variables and reoptimizes ALL continuous
variables in the original matrix. This improves upper bounds, never certifies
the original MILP optimum and never changes the original equations.
"""
from __future__ import annotations
import itertools
import time
import numpy as np
import highspy
from scipy.sparse import csr_matrix
from .solver import check_solution


class FixedPlanLP:
    def __init__(self,m,x):
        self.m=m
        m.solution=np.asarray(x,dtype=float).copy()
        residual=check_solution(m)
        if max(residual.values())>2e-6:raise ValueError(f'Invalid starting solution: {residual}')
        self.best=m.solution.copy();self.best_cost=m.value(m.objective)
        self.initial_cost=self.best_cost;self.tested=0;self.accepted=[];self.seen=set()
        self.integer=np.flatnonzero(m.integer).astype(np.int32)
        active=[r for r in m.rows if r.active];rr=[];cc=[];vv=[]
        for i,r in enumerate(active):
            for j,c in r.expression.terms.items():
                if c:rr.append(i);cc.append(j);vv.append(c)
        self.a=csr_matrix((vv,(rr,cc)),shape=(len(active),len(x)))
        self.row_lower=np.array([r.lower for r in active]);self.row_upper=np.array([r.upper for r in active])
        lower=np.asarray(m.lower).copy();upper=np.asarray(m.upper).copy()
        lower[self.integer]=np.round(x[self.integer]);upper[self.integer]=lower[self.integer]
        cost=np.zeros(len(x))
        for j,c in m.objective.terms.items():cost[j]=c
        lp=highspy.HighsLp();lp.num_col_=len(x);lp.num_row_=len(active)
        lp.col_cost_=cost;lp.col_lower_=lower;lp.col_upper_=upper;lp.offset_=m.objective.constant
        lp.row_lower_=self.row_lower;lp.row_upper_=self.row_upper
        lp.a_matrix_.format_=highspy.MatrixFormat.kRowwise
        lp.a_matrix_.start_=self.a.indptr;lp.a_matrix_.index_=self.a.indices;lp.a_matrix_.value_=self.a.data
        highspy.Highs.resetGlobalScheduler(True)
        self.h=highspy.Highs()
        for k,v in [('output_flag',False),('threads',1),('primal_feasibility_tolerance',1e-8)]:self.h.setOptionValue(k,v)
        self.h.passModel(lp)

    def trial(self,changes,label,deadline):
        remaining=deadline-time.perf_counter()
        if remaining<=0:return False
        values=np.round(self.best[self.integer]).astype(float)
        index_position={int(j):i for i,j in enumerate(self.integer)}
        for j,v in changes.items():values[index_position[int(j)]]=round(v)
        key=values.astype(np.int8).tobytes()
        if key in self.seen:return False
        self.seen.add(key)
        self.h.changeColsBounds(len(values),self.integer,values,values)
        self.h.setOptionValue('time_limit',self.h.getRunTime()+min(3.,remaining))
        self.h.run();self.tested+=1
        info=self.h.getInfo();sol=self.h.getSolution()
        if not sol.value_valid or info.primal_solution_status!=highspy.SolutionStatus.kSolutionStatusFeasible:return False
        if info.objective_function_value>=self.best_cost-1e-5:return False
        x=np.asarray(sol.col_value)
        activity=self.a@x
        residual=max(float(np.max(self.row_lower-activity,initial=0)),float(np.max(activity-self.row_upper,initial=0)))
        if residual>2e-6:return False
        self.best=x.copy();self.best_cost=float(info.objective_function_value)
        self.accepted.append(dict(candidate=label,objective=self.best_cost))
        print(f'JOINT {label}: {self.best_cost:.6f}',flush=True)
        return True


def search_joint_plans(m,p,s,x,seconds=120,strategies=('flat','swaps','timing','coordinate'),seed_solutions=()):
    evaluator=FixedPlanLP(m,x)
    start=time.perf_counter();deadline=start+max(0.,seconds)
    all_idx=evaluator.integer
    # Neighboring scenarios are only starting PLANS: all physical quantities are
    # repaired against the current model before a candidate can be accepted.
    for n,seed in enumerate(seed_solutions):
        if len(seed)==len(x):evaluator.trial(dict(zip(all_idx,np.round(seed[all_idx]))),f'neighbor {n}',min(deadline,start+seconds*.1))
    from .workflow import plan_values
    if 'flat' in strategies and np.ptp(p.demand.data)<1e-12:
        stage_end=deadline if len(strategies)==1 else min(deadline,start+seconds*.35)
        sites=[l for l in s.L if p.n_filter[l,5]]
        dates=[(7,8),(8,7),(7,7),(8,8),(6,8),(8,6),(6,7),(7,6),(6,6),(5,8),(8,5),(5,7),(7,5)]
        for ta,tb in dates:
            for a,b in itertools.combinations(sites,2):
                # Joint replacement at TWO sites can cross barriers that changing
                # only one site/one variable family at a time cannot cross.
                retire={a:ta,b:tb};new={a:(ta,5),b:(tb,5)}
                evaluator.trial(plan_values(m,retire=retire,new=new),f'new CCS {a}@{ta},{b}@{tb}',stage_end)
                for q in s.L:
                    if q in (a,b):continue
                    for tq in [3,4,5]:
                        evaluator.trial(plan_values(m,retire={**retire,q:tq},new=new),f'CCS pair; retire {q}@{tq}',stage_end)
                if time.perf_counter()>=stage_end:break
            if time.perf_counter()>=stage_end:break
    if 'swaps' in strategies:
        stage_end=deadline if len(strategies)==1 else min(deadline,time.perf_counter()+seconds*.10)
        for a,b in itertools.combinations(s.L,2):
            changes={}
            for family in ['y_o','y_e','y_r','y_n']:
                for key,v in m.variables[family].items():
                    if key[2] not in (a,b):continue
                    other=b if key[2]==a else a
                    donor=m.variables[family].get(key[:2]+(other,)+key[3:])
                    changes[v.index]=0. if donor is None else round(evaluator.best[donor.index])
            evaluator.trial(changes,f'swap sites {a},{b}',stage_end)
            if time.perf_counter()>=stage_end:break
    if 'timing' in strategies:
        stage_end=deadline if len(strategies)==1 else min(deadline,time.perf_counter()+seconds*.45)
        # A delayed retrofit at one site may be infeasible alone, while an earlier
        # retrofit at another site is more costly alone. Moving BOTH dates can
        # be feasible and cheaper: this escapes the coordinate-search barrier.
        for families in [('y_r',),('y_e',),('y_o','y_n'),('y_o','y_e','y_r','y_n')]:
            for a,b in itertools.combinations(s.L,2):
                for delta in [-1,1]:
                    changes={}
                    for family in families:
                        for key,v in m.variables[family].items():
                            if key[2] not in (a,b):continue
                            shift=delta if key[2]==a else -delta
                            t=(key[0]-1)*p.n_subperiods+key[1]
                            donor=min(p.n_periods*p.n_subperiods,t-shift)
                            if donor<1:
                                value=family=='y_o' or (family in ['y_r','y_n'] and key[3]==1)
                            else:
                                oldkey=((donor-1)//p.n_subperiods+1,(donor-1)%p.n_subperiods+1)+key[2:]
                                value=round(evaluator.best[m.variables[family][oldkey].index])
                            changes[v.index]=float(value)
                    evaluator.trial(changes,f'timing {families}:sites {a},{b},shift {delta}',stage_end)
                    if time.perf_counter()>=stage_end:break
                if time.perf_counter()>=stage_end:break
            if time.perf_counter()>=stage_end:break
    if 'coordinate' in strategies:
        for sweep in range(3):
            before=evaluator.best_cost
            for family in ['y_o','y_r','y_n','y_e']:
                for l in s.L:
                    entries=list((k,v.index) for k,v in m.variables[family].items() if k[2]==l)
                    if family in ['y_o','y_e']:
                        for onset in range(2 if family=='y_o' else 1,10):
                            changes={j:float((k[0]-1)*p.n_subperiods+k[1]<onset) if family=='y_o'
                                     else float((k[0]-1)*p.n_subperiods+k[1]>=onset) for k,j in entries}
                            evaluator.trial(changes,f'{family} site {l}@{onset}',deadline)
                    else:
                        for tech in sorted({k[3] for k,j in entries}):
                            for onset in ([9] if tech==1 else range(2,9)):
                                changes={j:float(k[3]==(tech if (k[0]-1)*p.n_subperiods+k[1]>=onset else 1)) for k,j in entries}
                                evaluator.trial(changes,f'{family} site {l}:tech{tech}@{onset}',deadline)
                    if time.perf_counter()>=deadline:break
                if time.perf_counter()>=deadline:break
            # Retirement and new construction at one site must also be movable
            # together, including independent replacement dates at other sites.
            for l in s.L:
                for tech in sorted({k[3] for k in m.variables['y_n'] if k[2]==l and k[3]>1}):
                    for onset in range(2,9):
                        changes={}
                        for family in ['y_o','y_n']:
                            for k,v in m.variables[family].items():
                                if k[2]!=l:continue
                                t=(k[0]-1)*p.n_subperiods+k[1]
                                changes[v.index]=float(t<onset) if family=='y_o' else float(k[3]==(tech if t>=onset else 1))
                        evaluator.trial(changes,f'replace site {l}:tech{tech}@{onset}',deadline)
                if time.perf_counter()>=deadline:break
            if evaluator.best_cost>=before-1e-5 or time.perf_counter()>=deadline:break
    m.solution=evaluator.best.copy()
    residual=check_solution(m)
    if max(residual.values())>2e-6:raise RuntimeError(f'Joint search produced invalid solution: {residual}')
    stats=dict(method='fixed-plan LP joint and coordinate neighborhoods; not global optimization',
        initial_objective=evaluator.initial_cost,final_objective=evaluator.best_cost,
        candidates=evaluator.tested,accepted=len(evaluator.accepted),improvements=evaluator.accepted,
        seconds=time.perf_counter()-start,global_lower_bound_updated=False,**residual)
    return evaluator.best.copy(),stats
