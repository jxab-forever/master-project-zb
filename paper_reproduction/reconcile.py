"""Reoptimize actual cross-target plans; no interpolation of curve values."""
import argparse
import json
import time
import numpy as np
from .workflow import ROOT, OUT, scenario_params, sets, build_softwarex, export_data
from .joint_search import FixedPlanLP
from .solver import check_solution


def repair_with_plans(m, incumbent, seeds, seconds=30):
    evaluator = FixedPlanLP(m, incumbent)
    started = time.perf_counter()
    deadline = started + seconds
    for n, seed in enumerate(seeds):
        if time.perf_counter() >= deadline:
            break
        seed = np.asarray(seed, dtype=float)
        if len(seed) != len(incumbent) or not np.isfinite(seed).all():
            continue
        # A complete solution from a stricter budget can already be feasible.
        # Check it directly so an LP timeout cannot discard that valid fallback.
        m.solution = seed
        if max(check_solution(m).values()) < 2e-6:
            cost = m.value(m.objective)
            if cost < evaluator.best_cost - 1e-5:
                evaluator.best = seed.copy()
                evaluator.best_cost = cost
                evaluator.accepted.append(dict(candidate=f'valid vector {n}', objective=cost))
        # Plans from a weaker budget must instead be repaired in this target's
        # complete original model before their quantities can be accepted.
        plan = dict(zip(evaluator.integer, np.round(seed[evaluator.integer])))
        evaluator.trial(plan, f'cross-target plan {n}', deadline)
    m.solution = evaluator.best.copy()
    residual = check_solution(m)
    if max(residual.values()) > 2e-6:
        raise RuntimeError(f'Invalid cross-target transfer: {residual}')
    stats = dict(method='complete-vector feasibility and fixed-plan LP cross-target repair',
        initial_objective=evaluator.initial_cost, final_objective=evaluator.best_cost,
        candidates=evaluator.tested, accepted=len(evaluator.accepted),
        improvements=evaluator.accepted, seconds=time.perf_counter()-started,
        global_lower_bound_updated=False, **residual)
    return evaluator.best.copy(), stats


def main(seconds=30, sweeps=2):
    for sweep in range(sweeps):
        order = range(9) if sweep % 2 == 0 else range(8, -1, -1)
        for i in order:
            name = f'target_{i:02d}'
            folder = OUT / 'refined' / name
            summary = json.loads((folder / 'summary.json').read_text())
            p = scenario_params(ROOT / 'inputs/softwarex/f_4.xlsx', .02, .025, 1, summary['target'])
            s = sets(p)
            m = build_softwarex(p, s)
            x = np.load(folder / 'solution.npy')
            others = sorted([j for j in range(9) if j != i], key=lambda j: abs(j-i))
            seeds = [np.load(OUT / 'refined' / f'target_{j:02d}' / 'solution.npy') for j in others]
            seeds.append(np.load(OUT / 'refined/nominal/solution.npy'))
            best, stats = repair_with_plans(m, x, seeds, seconds)
            old = summary['objective']
            value = m.value(m.objective)
            if value > old + 1e-6:
                raise RuntimeError('Reconciliation worsened a feasible incumbent')
            history = summary.get('refinement_passes', [summary['refinement']])
            summary.update(check_solution(m))
            summary.update(objective=value, primal_bound=value,
                mip_gap=max(0., (value-summary['dual_bound'])/abs(value)),
                refinement=stats, refinement_passes=history+[stats])
            export_data(m, p, s, folder, summary)
            print(f'RECONCILED {sweep+1} {name}: {old:.6f} -> {value:.6f}', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--seconds', type=float, default=30)
    parser.add_argument('--sweeps', type=int, default=2)
    args = parser.parse_args()
    main(args.seconds, args.sweeps)
