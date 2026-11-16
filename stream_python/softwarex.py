"""Solve and export the full SoftwareX industrial transition case."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from .softwarex_data import read_params, sets
from .softwarex_model import build_softwarex, gen_lp_window


def run_softwarex(input_path, output_dir, *, time_limit=300, mip_gap=0.05,
                  threads=4, warm_start_dir=None):
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    print('SoftwareX: reading Excel and building the complete MILP...',flush=True)
    started = time.perf_counter()
    p = read_params(input_path)
    s = sets(p)
    m = build_softwarex(p, s)
    build_seconds = time.perf_counter()-started
    print(f'SoftwareX: {len(m.column_names)} variables, model built in {build_seconds:.1f}s.',flush=True)
    starts = {}
    if warm_start_dir:
        for family, filename in [('y_o','yostate.csv'),('y_e','yestate.csv'),
                                  ('y_r','yrstate.csv'),('y_n','ynstate.csv')]:
            path = Path(warm_start_dir)/filename
            if not path.exists():
                continue
            for row in pd.read_csv(path).itertuples(index=False, name=None):
                key = tuple(int(x) for x in row[:-1])
                if key in m.variables[family]:
                    starts[m.variables[family][key].index] = float(row[-1])
    print(f'SoftwareX: solving with HiGHS, time limit {time_limit}s, {len(starts)} initial binary values.',flush=True)
    summary = m.solve(time_limit=time_limit, mip_gap=mip_gap, threads=threads,
                       log_file=output/'solver.log', model_file=output/'model.mps',
                       starts=starts, output=False)
    summary.update(build_seconds=build_seconds,
                   configuration=dict(time_limit=time_limit,mip_gap=mip_gap,threads=threads,warm_start=bool(starts)))
    return export_softwarex_results(m,p,s,input_path,output,summary)


def export_softwarex_results(m,p,s,input_path,output_dir,summary):
    output=Path(output_dir)
    value = m.value
    loan_window,pay_window=gen_lp_window(p,s)
    annual = []
    locations = []
    choices = []
    import openpyxl
    wb=openpyxl.load_workbook(input_path,read_only=True,data_only=True)
    labels={family:[str(row[0]) for row in list(wb[sheet].values)[1:]]
            for family,sheet in [('y_r','RF_label'),('y_n','NW_label')]}
    wb.close()
    cost_names = ['o_pay','n_pay','o_cfonm','o_cvonm','n_cfonm','n_cvonm',
                  't_ret_cost','o_u_cost','n_u_cost','o_ehf_cost','n_ehf_cost',
                  'o_ep1gcs_cost','n_ep1gcs_cost','o_x_in_cost','n_x_in_cost']
    costs = []
    for i in s.P:
        for j in s.P2:
            year = p.y0 + ((i-1)*p.n_subperiods+j-1)*p.yr_subperiod
            row = dict(period=i, subperiod=j, year=year, years=p.yr_subperiod,
                       demand_tonnes_per_year=p.demand[i,j]/p.sf_cap,
                       electricity_budget_MMBtu_per_year=p.elec_budget[i,j]/p.sf_elec,
                       capture_budget_tonnes_per_year=p.initial_co2*0.9/p.sf_em)
            for family in ['o_ep1ge','n_ep1ge','o_ups_e_mt_in','n_ups_e_mt_in',
                           'o_ep1gce','n_ep1gce','o_ep1gcs','n_ep1gcs']:
                row[family+'_tonnes_per_year'] = sum(value(m.variables[family][i,j,l]) for l in s.L)/p.sf_em
            row['co2_tonnes_per_year'] = sum(row[k+'_tonnes_per_year'] for k in
                          ['o_ep1ge','n_ep1ge','o_ups_e_mt_in','n_ups_e_mt_in'])
            row['existing_tonnes_per_year'] = sum(value(m.o_cp[i,j,l,p.key_node]) for l in s.L)/p.sf_cap
            row['new_tonnes_per_year'] = sum(value(m.n_cp[i,j,l,p.key_node]) for l in s.L)/p.sf_cap
            row['production_tonnes_per_year'] = row['existing_tonnes_per_year']+row['new_tonnes_per_year']
            row['electricity_MMBtu_per_year'] = sum(value(m.o_u[i,j,l,n]+m.n_u[i,j,l,n])
                                                    for l in s.L for n in s.Nd if p.nd_en_fltr[n])/p.sf_elec
            annual.append(row)
            cr = dict(year=year)
            for family in cost_names:
                cr[family+'_discounted_MUSD'] = sum(value(m.variables[family][i,j,l]) for l in s.L)*p.discount[i,j]/p.sf_cash
            cr['n_pay_discounted_MUSD']=p.discount[i,j]/p.sf_cash*sum(value(m.n_pay[i0,j0,l])
                for l in s.L for i0 in s.P for j0 in s.P2 if (i,j) in pay_window[i0,j0])
            cr['terminal_discounted_MUSD']=0.0
            if i==p.n_periods and j==p.n_subperiods:
                cr['terminal_discounted_MUSD']=1e-3*p.discount[i,j]/p.sf_cash*(
                    sum(value(m.o_loan_last[i,l,1]) for l in s.L)+
                    sum(value(m.n_loan[i,j,i0,j0,l]) for l in s.L for i0 in s.P for j0 in s.P2
                        if (i,j) in loan_window[i0,j0] and (i,j)!=(i0,j0)))
            costs.append(cr)
            for l in s.L:
                locations.append(dict(year=year, location=l,
                    existing_production_tonnes_per_year=value(m.o_cp[i,j,l,p.key_node])/p.sf_cap,
                    new_production_tonnes_per_year=value(m.n_cp[i,j,l,p.key_node])/p.sf_cap,
                    existing_online=value(m.y_o[i,j,l]), expansion=value(m.y_e[i,j,l])))
                for family in ['y_r','y_n']:
                    for k in (s.Kr if family=='y_r' else s.Kn):
                        if (i,j,l,k) in m.variables[family]:
                            choices.append(dict(year=year, location=l, form=family,
                                                technology=k,label=labels[family][k-1],
                                                activity=value(m.variables[family][i,j,l,k])))
    annual = pd.DataFrame(annual)
    total_co2 = float((annual.co2_tonnes_per_year*annual.years).sum())
    carbon_row = next(r for r in m.rows if r.active and r.name=='ag_co2_l_link_i_:()')
    carbon_budget = carbon_row.upper/p.sf_em
    summary.update(case='softwarex', solver='highspy / HiGHS',
        source_workbook=str(Path(input_path).resolve()),
        source_sha256=hashlib.sha256(Path(input_path).read_bytes()).hexdigest(),
        objective_MUSD=summary['objective']/p.sf_cash,
        total_co2_tonnes=total_co2, cumulative_co2_budget_tonnes=carbon_budget,
        minimum_demand_margin_tonnes_per_year=float((annual.production_tonnes_per_year-annual.demand_tonnes_per_year).min()),
        max_electricity_excess_MMBtu_per_year=float(np.maximum(annual.electricity_MMBtu_per_year-annual.electricity_budget_MMBtu_per_year,0).max()),
        max_capture_excess_tonnes_per_year=float(np.maximum(annual.o_ep1gce_tonnes_per_year+annual.o_ep1gcs_tonnes_per_year+
                         annual.n_ep1gce_tonnes_per_year+annual.n_ep1gcs_tonnes_per_year-annual.capture_budget_tonnes_per_year,0).max()),
        cost_accounting_error_MUSD=float(sum(sum(v for k,v in row.items() if k!='year') for row in costs)-summary['objective']/p.sf_cash),
        units=dict(cost='million USD',capacity='tonnes/year',emissions='tonnes CO2',
                   electricity='MMBtu/year',note='Workbook units sheet specifies ktonne/kMMBtu; exports undo sf_cap, sf_em and sf_elec.'))
    annual.to_csv(output/'annual_results.csv',index=False)
    pd.DataFrame(locations).to_csv(output/'location_results.csv',index=False)
    pd.DataFrame(choices).to_csv(output/'technology_choices.csv',index=False)
    pd.DataFrame(costs).to_csv(output/'discounted_costs.csv',index=False)
    pd.DataFrame(dict(variable=m.column_names, value=m.solution)).to_csv(output/'solution_variables.csv',index=False)
    (output/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')
    return summary
