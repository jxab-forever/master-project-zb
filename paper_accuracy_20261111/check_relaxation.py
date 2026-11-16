"""Measure the published Julia model's cost lower bound before enforcing integrality."""
import argparse
import json
import os
from pathlib import Path

import highspy
import numpy as np


def run(root, output):
    h = highspy.Highs()
    h.setOptionValue('output_flag', False)
    previous = Path.cwd()
    os.chdir(root)
    try:
        status = h.readModel('validation/revision_20261110/scenario_4_julia.mps')
    finally:
        os.chdir(previous)
    if status == highspy.HighsStatus.kError:
        raise ValueError('Unable to load independent Julia MPS')
    lp = h.getLp()
    indices = np.array([i for i, kind in enumerate(lp.integrality_) if int(kind)], dtype=np.int32)
    h.changeColsIntegrality(len(indices), indices, np.zeros(len(indices), dtype=np.uint8))
    h.setOptionValue('threads', 1)
    h.setOptionValue('time_limit', 90.)
    h.run()
    if h.getModelStatus() != highspy.HighsModelStatus.kOptimal:
        raise ValueError(h.modelStatusToString(h.getModelStatus()))
    values = np.asarray(h.getSolution().col_value)
    fractional = indices[np.abs(values[indices] - np.round(values[indices])) > 1e-6]
    groups = {}
    for i in fractional:
        group = lp.col_names_[i].split('[')[0]
        groups.setdefault(group, []).append(dict(variable=lp.col_names_[i], value=float(values[i])))
    summary = json.loads((root/'paper_reproduction/runs/nominal/summary.json').read_text(encoding='utf-8'))
    result = dict(source='scenario_4_julia.mps', relaxed_objective=h.getObjectiveValue(),
                  binary_variables=len(indices), fractional_decisions=len(fractional),
                  fractional_by_family={key: len(items) for key, items in groups.items()},
                  fractional_examples={key: items[:4] for key, items in groups.items()},
                  existing_incumbent=summary['objective'], existing_mip_lower_bound=summary['dual_bound'],
                  existing_mip_gap=summary['mip_gap'])
    output.write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('julia_relaxation.json'))
    args = parser.parse_args()
    run(args.root, args.output)
