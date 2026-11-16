"""Recombine facility schedules while retaining the published STREAM equations."""
import argparse
import itertools
import json
from pathlib import Path
import sys
import time

import highspy
import numpy as np


def schedule_bank(m, sites, folders, current):
    variables = [
        {(family, *key[:2], *key[3:]): var.index
         for family in ('y_o', 'y_e', 'y_r', 'y_n')
         for key, var in m.variables[family].items() if key[2] == site}
        for site in sites
    ]
    indices = [np.array(list(group.values()), dtype=np.int32) for group in variables]
    if sorted(np.concatenate(indices).tolist()) != np.flatnonzero(m.integer).tolist():
        raise ValueError('Facility groups must cover each binary variable exactly once')
    lower = np.asarray(m.lower)
    upper = np.asarray(m.upper)
    candidates = [{tuple(np.round(current[group]).astype(int)): 'incumbent'} for group in indices]
    for file in folders:
        donor = np.load(file)
        if donor.shape != current.shape:
            continue
        for site, group in enumerate(variables):
            for src in [site, *(i for i in range(len(indices)) if i != site)]:
                origin = variables[src]
                if any(key not in group and round(donor[index]) for key, index in origin.items()):
                    continue
                pattern = tuple(round(donor[origin[key]]) if key in origin else 0 for key in group)
                if np.any(pattern < lower[indices[site]]) or np.any(pattern > upper[indices[site]]):
                    continue
                candidates[site].setdefault(pattern, f'{file.parent.name}/site{src+1}')
    return indices, [np.asarray(list(bank), dtype=float) for bank in candidates]


def restricted_mip(evaluator, groups, bank, seconds):
    h = highspy.Highs()
    h.setOptionValue('output_flag', False)
    h.passModel(evaluator.h.getLp())
    m = evaluator.m
    for site, (indices, patterns) in enumerate(zip(groups, bank), start=1):
        if len(patterns) == 1:
            continue
        selectors = np.arange(h.getNumCol(), h.getNumCol() + len(patterns), dtype=np.int32)
        h.addVars(len(selectors), np.zeros(len(selectors)), np.ones(len(selectors)))
        h.changeColsIntegrality(len(selectors), selectors,
                               np.full(len(selectors), int(highspy.HighsVarType.kInteger), dtype=np.uint8))
        h.addRow(1, 1, len(selectors), selectors, np.ones(len(selectors)))
        for variable, states in zip(indices, patterns.T):
            if np.all(states == states[0]):
                continue
            h.changeColBounds(int(variable), float(m.lower[variable]), float(m.upper[variable]))
            used = selectors[states == 1]
            h.addRow(0, 0, len(used) + 1, np.r_[variable, used], np.r_[1.0, -np.ones(len(used))])
        for i, selector in enumerate(selectors):
            h.passColName(int(selector), f'site_{site}_schedule_{i}')
    binary = np.r_[evaluator.best, *[np.r_[1., np.zeros(len(p)-1)] for p in bank if len(p)>1]]
    h.setSolution(len(binary), np.arange(len(binary), dtype=np.int32), binary)
    h.setOptionValue('threads', 1)
    h.setOptionValue('time_limit', seconds)
    h.setOptionValue('mip_rel_gap', .001)
    h.setOptionValue('log_to_console', False)
    h.setOptionValue('output_flag', True)
    h.run()
    sol = h.getSolution()
    info = h.getInfo()
    if sol.value_valid and info.primal_solution_status == highspy.SolutionStatus.kSolutionStatusFeasible:
        x = np.asarray(sol.col_value[:len(evaluator.best)])
        evaluator.trial(dict(zip(evaluator.integer, x[evaluator.integer])),
                        'schedule MIP', time.perf_counter()+10)
    return dict(status=h.modelStatusToString(h.getModelStatus()), nodes=int(info.mip_node_count),
                restricted_lower_bound=float(info.mip_dual_bound), seconds=h.getRunTime())


def recombine(evaluator, indices, bank, seconds):
    deadline = time.perf_counter() + seconds
    pool = list(itertools.product(*(range(len(patterns)) for patterns in bank))) if np.prod(
        [len(patterns) for patterns in bank], dtype=np.int64) <= 3000 else None
    if pool is not None:
        for selected in pool:
            if time.perf_counter() >= deadline:
                break
            changes = {int(j): v for site, n in enumerate(selected) for j, v in zip(indices[site], bank[site][n])}
            evaluator.trial(changes, f'complete schedule combination {selected}', deadline)
    rng = np.random.default_rng(20261111)
    choices = [(site, n) for site, patterns in enumerate(bank) for n in range(len(patterns))]
    while time.perf_counter() < deadline:
        previous = evaluator.best_cost
        for k in rng.permutation(len(choices)):
            site, n = choices[k]
            evaluator.trial(dict(zip(indices[site], bank[site][n])), f'facility {site+1} schedule {n}', deadline)
        for _ in range(120):
            if time.perf_counter() >= deadline:
                break
            selected = rng.choice(len(indices), size=rng.choice([2, 2, 3, 4]), replace=False)
            patterns = [rng.integers(len(bank[s])) for s in selected]
            changes = {int(j): v for s, n in zip(selected, patterns) for j, v in zip(indices[s], bank[s][n])}
            evaluator.trial(changes, f'facility combination {list(zip(selected.tolist(), patterns))}', deadline)
        if pool is not None or (previous == evaluator.best_cost and len(evaluator.seen) > 3000):
            break


def timing_moves(evaluator, p, sites, seconds):
    deadline = time.perf_counter() + seconds
    periods = p.n_periods * p.n_subperiods
    variables = [{family: [(key, var.index) for key, var in evaluator.m.variables[family].items() if key[2] == site]
                  for family in ('y_o', 'y_e', 'y_r', 'y_n')} for site in sites]
    rng = np.random.default_rng(1485)

    def time_of(key):
        return (key[0] - 1) * p.n_subperiods + key[1]

    def get_events(group):
        result = []
        for family, entries in group.items():
            events = [(time_of(key), key[3] if len(key) == 4 else 0) for key, j in entries
                      if (round(evaluator.best[j]) == 0 if family == 'y_o' else round(evaluator.best[j]) == 1)
                      and (len(key) == 3 or key[3] != 1)]
            result.append(min(events, default=(periods+1, 1 if family in ('y_r', 'y_n') else 0)))
        return result

    def put_events(group, events):
        values = {}
        for (family, entries), (date, technology) in zip(group.items(), events):
            for key, j in entries:
                active = time_of(key) >= date
                values[j] = int(not active if family == 'y_o' else active if family == 'y_e'
                                else key[3] == (technology if active else 1))
        return values

    def variants(group):
        original = get_events(group)
        if any(round(evaluator.best[j]) != v for j, v in put_events(group, original).items()):
            return []
        off, expansion, retrofit, replacement = original
        result = set()
        for i, (family, entries) in enumerate(group.items()):
            date, technology = original[i]
            for shift in (-2, -1, 1, 2):
                new_date = max(1 if family == 'y_e' else 2, min(periods+1, date+shift))
                changed = original.copy()
                changed[i] = (new_date, technology)
                result.add(tuple(changed))
            if family in ('y_r', 'y_n'):
                for tech in sorted({key[3] for key, _ in entries}):
                    for date_variant in ({date, max(2, date-1), periods} if tech != 1 else {periods+1}):
                        changed = original.copy()
                        changed[i] = (date_variant, tech)
                        result.add(tuple(changed))
        for delta in (-2, -1, 1, 2):
            new_off = max(2, min(periods+1, off[0] + delta))
            new_build = max(new_off, min(periods+1, replacement[0] + delta))
            new_retrofit = (new_off, retrofit[1]) if retrofit[0] == off[0] else retrofit
            result.add(((new_off, off[1]), expansion, new_retrofit, (new_build, replacement[1])))
        return [put_events(group, event) for event in sorted(result)
                if event[0][0] <= event[3][0]
                and (event[2][1] == 1 or event[2][0] <= event[0][0])
                and event != tuple(original)]

    while time.perf_counter() < deadline:
        pool = [variants(group) for group in variables]
        choices = [(site, option) for site, options in enumerate(pool) for option in range(len(options))]
        if not choices:
            break
        for n in rng.permutation(len(choices)):
            if time.perf_counter() >= deadline:
                break
            site, option = choices[n]
            evaluator.trial(pool[site][option], f'facility {site+1} event change {option}', deadline)
        for _ in range(60):
            if time.perf_counter() >= deadline:
                break
            a, b = rng.choice(len(choices), size=2, replace=False)
            sa, ia = choices[a]
            sb, ib = choices[b]
            if sa == sb:
                continue
            evaluator.trial({**pool[sa][ia], **pool[sb][ib]},
                            f'linked event changes at facilities {sa+1},{sb+1}', deadline)


def flat_capacity_seeds(evaluator, p, sites, seconds):
    from paper_reproduction.workflow import plan_values

    if np.ptp(p.demand.data) > 1e-10:
        raise ValueError('Flat capacity starts require constant production demand')
    deadline = time.perf_counter() + seconds
    # The capacity drop in Figure 6(e–h) is near the smallest site's capacity.
    # Use it to order starts, then solve and check against the unchanged model.
    ordered = sorted(sites, key=lambda site: abs(p.c0[site] / p.sf_cap - 69.4))
    ccs_sites = [site for site in sites if p.n_filter[site, 5]]
    for first in ordered[:2]:
        for start in (4, 3, 5, 6):
            early_off = {var.index: int((key[0]-1)*p.n_subperiods+key[1] < start)
                         for key, var in evaluator.m.variables['y_o'].items() if key[2] == first}
            evaluator.trial(early_off, f'flat demand incumbent with facility {first} retiring at {start}', deadline)
    for first in ordered[:2]:
        for a, b in itertools.combinations(ccs_sites, 2):
            for a_time, b_time in ((8, 8), (7, 8), (8, 7)):
                build = {a: (a_time, 5), b: (b_time, 5)}
                retirement = {a: a_time, b: b_time, first: 4}
                other_sites = [site for site in sites if site not in retirement]
                for other in [None, *other_sites]:
                    retrofit = {other: (8, 2)} if other else {}
                    for start in [4, 3, 5]:
                        if time.perf_counter() >= deadline:
                            return
                        evaluator.trial(plan_values(evaluator.m, retire={**retirement, first: start},
                                                    new=build, retro=retrofit),
                                        f'flat demand retire {first}@{start}; CCS {a}@{a_time},{b}@{b_time}; retrofit {other}',
                                        deadline)


def solve(args):
    sys.path.insert(0, str(args.root))
    from paper_reproduction.workflow import OUT, FACTORS, build_softwarex, export_data, scenario_params, sets
    from paper_reproduction.joint_search import FixedPlanLP
    from paper_reproduction.solver import check_solution
    from paper_reproduction.revision_audit import validate_case

    class PresolvedPlanLP(FixedPlanLP):
        def trial(self, changes, label, deadline):
            self.h.clearSolver()
            return super().trial(changes, label, deadline)

    revision = OUT / 'revision_20261110'
    sources = [OUT / 'runs', OUT / 'refined', OUT / 'initial_milp', revision / 'first_pass',
               revision / 'before_runs', args.output / 'runs']
    output = args.output / 'runs' / args.case
    output.mkdir(parents=True, exist_ok=True)
    source = output if (output / 'solution.npy').is_file() else OUT / 'runs' / args.case
    summary = json.loads((source / 'summary.json').read_text(encoding='utf-8'))
    target = summary['target']
    factors = (summary['demand_growth'], summary['electricity_growth'], summary['ccs_cost_factor'])
    p = scenario_params(args.root / 'inputs/softwarex/f_4.xlsx', *factors, target)
    s = sets(p)
    model = build_softwarex(p, s)
    if summary['bau']:
        model.remove_group('ag_co2_l_link_i_')
    evaluator_type = PresolvedPlanLP if args.cold_presolve else FixedPlanLP
    evaluator = evaluator_type(model, np.load(source / 'solution.npy'))
    files = [path for directory in sources for path in sorted(directory.glob('*/solution.npy'))]
    files.extend(sorted(revision.glob('*/solution.npy')))
    indices, bank = schedule_bank(model, s.L, files, evaluator.best)
    print(args.case, 'initial', evaluator.initial_cost, 'facility pools', [len(p) for p in bank], flush=True)
    mip = restricted_mip(evaluator, indices, bank, args.mip_seconds) if args.mip_seconds else None
    seed_deadline = time.perf_counter() + 3*len(args.seed)
    for seed in args.seed:
        values = np.load(seed)
        if values.shape != evaluator.best.shape:
            raise ValueError(f'Unexpected solution shape in {seed}')
        evaluator.trial(dict(zip(evaluator.integer, values[evaluator.integer])),
                        f'complete plan from {seed.parent.name}', seed_deadline)
    if args.flat_seconds:
        flat_capacity_seeds(evaluator, p, s.L, args.flat_seconds)
    recombine(evaluator, indices, bank, args.lp_seconds)
    timing_moves(evaluator, p, s.L, args.timing_seconds)
    model.solution = evaluator.best
    residuals = check_solution(model)
    if max(residuals.values()) > 2e-6:
        raise ValueError(residuals)
    history = summary.get('schedule_searches', [])
    history.append(dict(initial_objective=evaluator.initial_cost, final_objective=evaluator.best_cost,
                        plans_by_facility=[len(p) for p in bank], candidates=evaluator.tested,
                        improvements=evaluator.accepted, restricted_mip=mip,
                        cold_presolve=args.cold_presolve, explicit_seeds=[str(path) for path in args.seed]))
    summary.update(**residuals, objective=evaluator.best_cost, primal_bound=evaluator.best_cost,
                   mip_gap=max(0., (evaluator.best_cost-summary['dual_bound'])/abs(evaluator.best_cost)),
                   status='Feasible incumbent; not globally certified', schedule_searches=history)
    export_data(model, p, s, output, summary)
    reference = args.root / f'validation/revision_20261110/scenario_{FACTORS.index(factors)+1}_julia.mps'
    audit = validate_case(output, reference, target, summary['bau'], factors)
    (output / 'independent_check.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
    if not audit['passed']:
        raise ValueError(audit)
    print(args.case, 'verified', evaluator.initial_cost, '->', evaluator.best_cost, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case')
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--mip-seconds', type=float, default=60)
    parser.add_argument('--lp-seconds', type=float, default=45)
    parser.add_argument('--timing-seconds', type=float, default=0)
    parser.add_argument('--flat-seconds', type=float, default=0)
    parser.add_argument('--cold-presolve', action='store_true')
    parser.add_argument('--seed', type=Path, action='append', default=[])
    solve(parser.parse_args())
