"""Development-time parity audit. Never needed to run either Python case."""
from __future__ import annotations
import json
from pathlib import Path
import re
from validation.reference_mps import load_reference_mps
import numpy as np
import highspy
from scipy.sparse import csc_matrix, csr_matrix
from stream_python.softwarex_data import read_params, sets
from stream_python.softwarex_model import build_softwarex


def canonical(name):
    if ':' in name:
        group, key = name.split(':', 1)
        name = group + '[' + key + ']'
    name = re.sub(r'[\s()]', '', name)
    name = re.sub(r'(?<=[\[,])_(?=\d)', '', name)
    return name.replace(',]', ']').replace('[]','')


def compare(root):
    root = Path(root)
    p = read_params(root/'inputs/softwarex/f_4.xlsx')
    m = build_softwarex(p, sets(p))
    h = highspy.Highs()
    h.setOptionValue('output_flag', False)
    reference=root/'validation/revision_20261110/scenario_4_julia.mps'
    load_reference_mps(h,reference)
    lp = h.getLp()
    julia_col = {canonical(n): i for i,n in enumerate(lp.col_names_)}
    julia_row = {canonical(n): i for i,n in enumerate(lp.row_names_)}
    py_col = {canonical(n):i for i,n in enumerate(m.column_names)}
    active = [r for r in m.rows if r.active]
    py_row = {canonical(r.name):i for i,r in enumerate(active)}
    result = dict(julia_variables=lp.num_col_,python_variables=len(py_col),
                  julia_constraints=lp.num_row_,python_constraints=len(py_row))
    result['missing_python_variables'] = sorted(set(julia_col)-set(py_col))[:30]
    result['extra_python_variables'] = sorted(set(py_col)-set(julia_col))[:30]
    result['missing_python_constraints'] = sorted(set(julia_row)-set(py_row))[:30]
    result['extra_python_constraints'] = sorted(set(py_row)-set(julia_row))[:30]
    if set(julia_col)!=set(py_col) or set(julia_row)!=set(py_row):
        result['passed'] = False
    else:
        cols = np.array([julia_col[canonical(n)] for n in m.column_names])
        rows = np.array([julia_row[canonical(r.name)] for r in active])
        matrix = csc_matrix((lp.a_matrix_.value_,lp.a_matrix_.index_,lp.a_matrix_.start_),
                            shape=(lp.num_row_,lp.num_col_)).tocsr()[rows,:][:,cols]
        rr, cc, vv = [],[],[]
        for i,r in enumerate(active):
            for c,v in r.expression.terms.items():
                if v: rr.append(i); cc.append(c); vv.append(v)
        python = csr_matrix((vv,(rr,cc)),shape=matrix.shape)
        differences = (matrix-python).tocoo()
        result['max_coefficient_difference'] = float(np.max(np.abs(differences.data),initial=0))
        mismatched = np.flatnonzero(np.abs(differences.data)>1e-8)
        result['coefficient_mismatch_count'] = int(len(mismatched))
        result['coefficient_examples'] = [dict(row=active[int(differences.row[k])].name,
            column=m.column_names[int(differences.col[k])],difference=float(differences.data[k]))
            for k in mismatched[:20]]
        for label, julia_values, python_values in [
            ('column_lower', np.asarray(lp.col_lower_)[cols],np.asarray(m.lower)),
            ('column_upper', np.asarray(lp.col_upper_)[cols],np.asarray(m.upper)),
            ('row_lower', np.asarray(lp.row_lower_)[rows],np.array([r.lower for r in active])),
            ('row_upper', np.asarray(lp.row_upper_)[rows],np.array([r.upper for r in active]))]:
            equal = np.isclose(julia_values,python_values,rtol=1e-10,atol=1e-8)
            result[label+'_mismatch_count'] = int(np.count_nonzero(~equal))
            result[label+'_examples'] = [dict(name=(active[k].name if label.startswith('row') else m.column_names[k]),
                julia=float(julia_values[k]),python=float(python_values[k])) for k in np.flatnonzero(~equal)[:20]]
        cost = np.zeros(len(cols))
        for c,v in m.objective.terms.items(): cost[c] = v
        result['max_objective_difference'] = float(np.max(np.abs(np.asarray(lp.col_cost_)[cols]-cost),initial=0))
        result['objective_offset_difference'] = float(lp.offset_-m.objective.constant)
        integers = np.array([int(v) for v in lp.integrality_])[cols] != 0
        result['integrality_mismatch_count'] = int(np.count_nonzero(integers!=m.integer))
        result['passed'] = not any(v for k,v in result.items() if k.endswith('_mismatch_count')) and result['max_objective_difference']<1e-8
    target=root/'validation/softwarex_parity.json'
    target.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
    return result


if __name__=='__main__':
    compare(Path(__file__).parents[1])
