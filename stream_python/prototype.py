"""Excel-driven Python port of the original continuous STRE3AM prototype.

The LP is rebuilt from the workbook, not read from an exported Julia model.
Purely definitional JuMP variables (generation, heat, cost, and effective
capacity) are substituted algebraically. Asset stocks, retirement flows,
retrofit transitions, and construction allocations retain the original timing.

Port of Argonne's ``mat_struct.jl``, ``coef_custom.jl``, ``model.jl``, and
``modKern.jl``. Copyright (C) 2023, UChicago Argonne, LLC. All Rights Reserved.
Software Name: STRE3AM: Strategic Technology Roadmapping and Energy,
Environmental, and Economic Analysis Model. By: Argonne National Laboratory.
BSD-3 OPEN SOURCE LICENSE.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from openpyxl import load_workbook
from scipy.optimize import linprog
from scipy import __version__ as scipy_version
from scipy.sparse import coo_matrix, csr_matrix


SHORT_NAMES = ("Coal", "GasCT", "GasCC", "Oil", "Biomass", "Nuclear",
               "Hydro", "Wind", "SolarPV", "SolarTH", "Geo")
COST_CATEGORIES = ("capital", "fixed_om", "variable_om", "fuel", "retirement")


@dataclass(frozen=True)
class Form:
    delay: int
    service_increase: float
    capital: float
    capital_base: int
    variable: float
    variable_base: int
    fixed: float
    fixed_base: int
    heat: float
    heat_base: int
    emissions: float
    emissions_base: int
    fuel: float
    fuel_base: int
    use_base_heat: bool
    ccs_efficiency_penalty: bool
    efficiency_points: float


@dataclass
class PrototypeData:
    path: Path
    T: int
    I: int
    names: list[str]
    initial: np.ndarray
    demand: np.ndarray
    capacity_factor: np.ndarray
    old_heat: np.ndarray
    new_heat: np.ndarray
    capital: np.ndarray
    fixed: np.ndarray
    variable: np.ndarray
    sale: np.ndarray
    fuel: np.ndarray
    decommission: np.ndarray
    service_life: np.ndarray
    carbon: np.ndarray
    Kz: np.ndarray
    Kx: np.ndarray
    fuel_based: np.ndarray
    co2_based: np.ndarray
    baseload: np.ndarray
    N: np.ndarray
    discount_rate: float
    heat_increase: float
    loan_period: int
    retrofit: dict[tuple[int, int], Form]
    new: dict[tuple[int, int], Form]
    generation_scale: float
    heat_to_efficiency: float

    @classmethod
    def from_excel(cls, path: str | Path) -> "PrototypeData":
        path = Path(path).resolve()
        book = load_workbook(path, data_only=True, read_only=True)
        try:
            ref = book["reference"]

            def read(sheet: str, row: int, factor: bool = False) -> np.ndarray:
                address = str(ref.cell(row, 2).value)
                cells = book[sheet][address]
                if not isinstance(cells, tuple):
                    cells = ((cells,),)
                elif cells and not isinstance(cells[0], tuple):
                    cells = (cells,)
                values = np.array([[c.value for c in r] for r in cells], dtype=float)
                if not np.all(np.isfinite(values)):
                    raise ValueError(f"Missing/non-numeric cached values in {sheet}!{address}")
                return values * (float(ref.cell(row, 3).value) if factor else 1.0)

            inv = book["invrAttr"]
            I, T = int(inv["B1"].value), int(inv["B2"].value)
            Kz = read("invrAttr", 22).astype(int).ravel(order="F")
            Kx = read("invrAttr", 23).astype(int).ravel(order="F")

            def forms(kinds: np.ndarray, ref_row: int) -> dict[tuple[int, int], Form]:
                cells = inv[str(ref.cell(ref_row, 2).value)]
                table = [[cell.value for cell in row] for row in cells]
                result: dict[tuple[int, int], Form] = {}
                offset = 0
                for i, count in enumerate(kinds):
                    for k in range(count):
                        v = table[offset + 1 + k]
                        result[i, k] = Form(
                            max(int(v[0]), 0), max(float(v[1]), 0),
                            *[float(v[j]) if j % 2 == 0 else math.floor(v[j])
                              for j in range(2, 14)],
                            bool(v[14]), bool(v[15]), float(v[16] or 0))
                    offset += int(count) + 1
                return result

            return cls(
                path, T, I, [str(inv.cell(11 + i, 1).value) for i in range(I)],
                read("timeAttr", 2, True), read("timeAttr", 3, True).ravel(order="F"),
                read("timeAttr", 4), read("timeAttr", 5, True),
                read("timeAttr", 6, True), read("costAttr", 8, True),
                read("costAttr", 9, True), read("costAttr", 10, True),
                read("costAttr", 11, True).ravel(order="F"),
                read("costAttr", 12, True), read("costAttr", 13, True),
                read("invrAttr", 17).astype(int).ravel(order="F"),
                read("invrAttr", 18, True).ravel(order="F"), Kz, Kx,
                read("invrAttr", 24).astype(bool).ravel(order="F"),
                read("invrAttr", 25).astype(bool).ravel(order="F"),
                read("invrAttr", 26).astype(bool).ravel(order="F"),
                read("invrAttr", 27).astype(int).ravel(order="F"),
                float(read("invrAttr", 19).item()),
                float(read("invrAttr", 20).item()),
                int(read("invrAttr", 21).item()), forms(Kz, 28), forms(Kx, 29),
                float(ref["B31"].value), float(ref["B32"].value))
        finally:
            book.close()

    def discount(self, t: int) -> float:
        return (1.0 + self.discount_rate) ** (-t)

    def retirement_cost(self, i: int, k: int, t: int, j: int, tag: str = "") -> float:
        age = t - j if tag == "X" else t + j
        base_age = min(j, self.capital.shape[1] - 1) if tag == "X" else 0
        multiplier = 1.0
        if tag == "R":
            multiplier = self.retrofit[i, k].capital
            if multiplier < 1:
                multiplier += 1
        elif tag == "X":
            multiplier = self.new[i, k].capital
        return (max(self.loan_period - age, 0) / self.loan_period * multiplier
                * self.capital[i, base_age] * self.discount(t))

    def lost_sales(self, i: int, t: int, j: int, tag: str = "") -> float:
        age = t - j if tag == "X" else t + j
        return max(self.service_life[i] - age, 0) * self.sale[t] * self.discount(t)

    def heat_rate(self, form: str, i: int, k: int, j: int, t: int) -> float:
        if form == "existing":
            return self.old_heat[i, j] * (1 + self.heat_increase) ** t
        f = (self.retrofit if form == "retrofit" else self.new)[i, k]
        base = f.heat_base if f.heat_base >= 0 else i
        hr = self.new_heat[base, min(j, self.new_heat.shape[1] - 1)]
        multiplier = f.heat
        if form == "retrofit":
            if f.use_base_heat:
                hr = self.old_heat[base, j]
                if hr <= 1e-6:
                    positive = self.old_heat[base][self.old_heat[base] > 0]
                    hr = float(positive.mean()) if len(positive) else 0.0
            if f.ccs_efficiency_penalty:
                multiplier = 1.0
                efficiency = np.inf if hr == 0 else 1 / (hr * self.heat_to_efficiency)
                efficiency -= f.efficiency_points / 100
                if efficiency > 0:
                    hr = 1 / efficiency / self.heat_to_efficiency
        return multiplier * hr * (1 + self.heat_increase) ** t

    def derating(self, i: int, k: int, j: int) -> float:
        f = self.retrofit[i, k]
        if not f.ccs_efficiency_penalty:
            return 1.0
        hr = self.old_heat[i, j]
        efficiency = np.inf if hr == 0 else 1 / (hr * self.heat_to_efficiency)
        if efficiency <= 0:
            return 1.0
        multiplier = 1 - f.efficiency_points / 100 / efficiency
        return multiplier if 0 <= multiplier <= 1 else 1e-8


@dataclass
class PrototypeLP:
    data: PrototypeData
    variables: dict[str, dict[tuple[int, ...], int]]
    objective: np.ndarray
    npv: np.ndarray
    latent: np.ndarray
    terminal: np.ndarray
    A_eq: csr_matrix
    b_eq: np.ndarray
    A_ub: csr_matrix
    b_ub: np.ndarray
    costs: csr_matrix
    emissions: csr_matrix
    # One row per (year, technology, form); the columns are LP variables.
    capacity: csr_matrix
    generation: csr_matrix
    asset_keys: list[tuple[int, int, str]]
    build_seconds: float


def build_prototype(data: PrototypeData | str | Path, *, latent_factor: float = 1e-3,
                    no_retrofit: bool = False, no_delay_z: bool = False,
                    no_delay_x: bool = False) -> PrototypeLP:
    """Build the mathematical model, eliminating only definitional auxiliaries.

    Original quirks are intentional: year 0 has no demand lower bound; retrofit
    fixed/variable/capital costs ignore their base-fuel selectors; retrofit
    latent life uses the child index; retrofit terminal lost sales omit the
    capacity factor; all effective capacities use the post-retirement stock.
    """
    started = time.perf_counter()
    d = data if isinstance(data, PrototypeData) else PrototypeData.from_excel(data)
    T, I = d.T, d.I
    kz = np.zeros(I, dtype=int) if no_retrofit else d.Kz
    dz = {(i, k): 0 if no_delay_z else d.retrofit[i, k].delay
          for i in range(I) for k in range(kz[i])}
    dx = {(i, k): 0 if no_delay_x else d.new[i, k].delay
          for i in range(I) for k in range(d.Kx[i])}
    groups: dict[str, dict[tuple[int, ...], int]] = {
        name: {} for name in ("w", "uw", "z", "uz", "zTrans", "x", "ux", "xAlloc")}
    count = 0

    def var(group: str, key: tuple[int, ...]) -> int:
        nonlocal count
        groups[group][key] = count
        count += 1
        return count - 1

    # The original's unused time-T allocations, retirements and new cohorts
    # have no constraints/objective contribution and can be omitted safely.
    for i in range(I):
        for j in range(d.N[i]):
            for t in range(T + 1):
                var("w", (t, i, j))
                if t < T:
                    var("uw", (t, i, j))
            for k in range(kz[i]):
                for t in range(T):
                    var("zTrans", (t, i, k, j))
                for t in range(dz[i, k], T + 1):
                    var("z", (t, i, k, j))
                    if t < T:
                        var("uz", (t, i, k, j))
        for k in range(d.Kx[i]):
            for t in range(T):
                var("xAlloc", (t, i, k))
            for j in range(T):
                for t in range(j + dx[i, k], T + 1):
                    var("x", (t, i, k, j))
                    if t < T:
                        var("ux", (t, i, k, j))

    eqr: list[int] = []
    eqc: list[int] = []
    eqv: list[float] = []
    eqb: list[float] = []

    def equation(terms: list[tuple[int, float]], rhs: float = 0) -> None:
        row = len(eqb)
        for index, coefficient in terms:
            if coefficient:
                eqr.append(row)
                eqc.append(index)
                eqv.append(coefficient)
        eqb.append(rhs)

    w, uw, z, uz, zt, x, ux, xa = (groups[g] for g in groups)
    for i in range(I):
        for j in range(d.N[i]):
            equation([(w[0, i, j], 1)], d.initial[i, j])
            for t in range(T):
                equation([(w[t + 1, i, j], 1), (w[t, i, j], -1),
                          (uw[t, i, j], 1)] +
                         [(zt[t, i, k, j], 1) for k in range(kz[i])])
            for k in range(kz[i]):
                delay = dz[i, k]
                if delay > T:
                    raise ValueError("Original model requires retrofit delay <= horizon")
                equation([(z[delay, i, k, j], 1)])
                for t in range(delay, T):
                    equation([(z[t + 1, i, k, j], 1), (z[t, i, k, j], -1),
                              (uz[t, i, k, j], 1), (zt[t - delay, i, k, j], -1)])
        for k in range(d.Kx[i]):
            for j in range(T):
                delay = dx[i, k]
                if j + delay < T:
                    equation([(x[j + delay, i, k, j], 1), (xa[j, i, k], -1)])
                for t in range(j + delay, T):
                    equation([(x[t + 1, i, k, j], 1), (x[t, i, k, j], -1),
                              (ux[t, i, k, j], 1)])

    costr: list[int] = []
    costc: list[int] = []
    costv: list[float] = []
    npv = np.zeros(count)
    latent = np.zeros(count)
    terminal = np.zeros(count)

    def cost(t: int, category: str, index: int, coefficient: float) -> None:
        if coefficient:
            npv[index] += coefficient
            costr.append(t * len(COST_CATEGORIES) + COST_CATEGORIES.index(category))
            costc.append(index)
            costv.append(coefficient)

    # Capital costs are charged on the allocation/transition date.
    for (t, i, k), index in xa.items():
        cost(t, "capital", index, d.new[i, k].capital * d.capital[i, t] * d.discount(t))
    for (t, i, k, j), index in zt.items():
        cost(t, "capital", index,
             d.retrofit[i, k].capital * d.capital[i, t] * d.discount(t))

    capr: list[int] = []
    capc: list[int] = []
    capv: list[float] = []
    genr: list[int] = []
    genc: list[int] = []
    genv: list[float] = []
    emr: list[int] = []
    emc: list[int] = []
    emv: list[float] = []
    asset_keys = [(t, i, form) for t in range(T) for i in range(I)
                  for form in ("existing", "retrofit", "new")]
    asset_rows = {key: row for row, key in enumerate(asset_keys)}

    def active(form: str, t: int, i: int, k: int, j: int, index: int) -> None:
        row = asset_rows[t, i, form]
        f = None if form == "existing" else (d.retrofit if form == "retrofit" else d.new)[i, k]
        cf = d.capacity_factor[i, j if form == "new" else 0]
        gc = d.generation_scale * 8760 * cf
        if form == "retrofit":
            gc *= d.derating(i, k, j)
        capr.append(row)
        capc.append(index)
        capv.append(1.0)
        genr.append(row)
        genc.append(index)
        genv.append(gc)
        discount = d.discount(t)
        cost(t, "fixed_om", index, d.fixed[i, t] * discount * (1 if f is None else f.fixed))
        cost(t, "variable_om", index,
             d.variable[i, t] * discount * gc * (1 if f is None else f.variable))
        if d.fuel_based[i]:
            heat = gc * d.heat_rate(form, i, k, j, t)
            fuel_base = i if f is None or f.fuel_base < 0 else f.fuel_base
            cost(t, "fuel", index, heat * d.fuel[fuel_base, t] * discount
                 * (1 if f is None else f.fuel))
            if d.co2_based[i]:
                carbon_base = i if f is None or f.emissions_base < 0 else f.emissions_base
                emr.append(t)
                emc.append(index)
                emv.append(heat * d.carbon[carbon_base] * (1 if f is None else f.emissions))
        if t == T - 1:
            if form == "new":
                tc = d.retirement_cost(i, k, t, j, "X") + (
                    d.generation_scale * 8760 * cf * d.lost_sales(i, t, j, "X"))
            elif form == "retrofit":
                # Source termCzE has neither tag="R" nor the capacity factor.
                tc = d.retirement_cost(i, k, t, j) + (
                    d.generation_scale * 8760 * d.lost_sales(i, t, j))
            else:
                tc = d.retirement_cost(i, i, t, j) + (
                    d.generation_scale * 8760 * cf * d.lost_sales(i, t, j))
            terminal[index] += tc

    for i in range(I):
        for j in range(d.N[i]):
            for t in range(T):
                active("existing", t, i, 0, j, w[t + 1, i, j])
                cost(t, "retirement", uw[t, i, j], d.retirement_cost(i, i, t, j) +
                     d.generation_scale * 8760 * d.capacity_factor[i, 0] * d.lost_sales(i, t, j))
                if t + j >= d.service_life[i] and not d.baseload[i]:
                    latent[w[t, i, j]] += d.retirement_cost(0, 0, 1, 0) * math.exp(
                        (j + t - d.service_life[i]) / d.service_life[i])
            for k in range(kz[i]):
                for t in range(dz[i, k], T):
                    active("retrofit", t, i, k, j, z[t + 1, i, k, j])
                    cost(t, "retirement", uz[t, i, k, j], d.retirement_cost(i, k, t, j, "R") +
                         d.generation_scale * 8760 * d.capacity_factor[i, 0] * d.lost_sales(i, t, j))
                    if t + j >= d.service_life[i] and not d.baseload[i]:
                        # Source rLatentZ uses servLife[k+1], not servLife[i+1].
                        life = d.service_life[k]
                        latent[z[t, i, k, j]] += d.retirement_cost(0, 0, 1, 0) * math.exp(
                            (j + t - life) / life)
        for k in range(d.Kx[i]):
            for j in range(T):
                for t in range(j + dx[i, k], T):
                    active("new", t, i, k, j, x[t + 1, i, k, j])
                    cost(t, "retirement", ux[t, i, k, j], d.retirement_cost(i, k, t, j, "X") +
                         d.generation_scale * 8760 * d.capacity_factor[i, j] * d.lost_sales(i, t, j, "X"))
                    if t - j >= d.service_life[i] and not d.baseload[i]:
                        life = math.floor(d.service_life[i] * (1 + d.new[i, k].service_increase))
                        latent[x[t, i, k, j]] += d.retirement_cost(0, 0, 1, 0) * math.exp(
                            (t - j - life) / life)

    shape = (len(asset_keys), count)
    capacity = coo_matrix((capv, (capr, capc)), shape=shape).tocsr()
    generation = coo_matrix((genv, (genr, genc)), shape=shape).tocsr()
    emissions = coo_matrix((emv, (emr, emc)), shape=(T, count)).tocsr()
    # Demand starts at t=1 in the original formulation.
    aggregate = coo_matrix((np.ones(len(asset_keys)),
                           ([t for t, _, _ in asset_keys], range(len(asset_keys)))),
                          shape=(T, len(asset_keys))).tocsr()
    generation_by_year = aggregate @ generation
    A_ub = -generation_by_year[1:]
    b_ub = -d.demand[1:T]
    return PrototypeLP(
        d, groups, (npv + latent_factor * latent + 1e-6 * terminal) / 1e3,
        npv, latent, terminal,
        coo_matrix((eqv, (eqr, eqc)), shape=(len(eqb), count)).tocsr(),
        np.array(eqb), A_ub, b_ub,
        coo_matrix((costv, (costr, costc)), shape=(T * len(COST_CATEGORIES), count)).tocsr(),
        emissions, capacity, generation, asset_keys, time.perf_counter() - started)


def _write_csv(path: Path, header: list[str], rows: Any) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def run_prototype(input_path: str | Path, output_dir: str | Path,
                  *, latent_factor: float = 1e-3, no_retrofit: bool = False,
                  no_delay_z: bool = False, no_delay_x: bool = False) -> dict[str, Any]:
    """Solve the workbook using SciPy/HiGHS and export consistent physical data.

    Units after the workbook's reference factors: GW, TWh/year, Mt CO2/year,
    and discounted million USD. The objective divides augmented cost by 1000.
    Capacity exports are *effective operating* capacity, including age zero.
    """
    started = time.perf_counter()
    model = build_prototype(input_path, latent_factor=latent_factor,
                            no_retrofit=no_retrofit, no_delay_z=no_delay_z,
                            no_delay_x=no_delay_x)
    solve_start = time.perf_counter()
    result = linprog(model.objective, A_ub=model.A_ub, b_ub=model.b_ub,
                     A_eq=model.A_eq, b_eq=model.b_eq, bounds=(0, None),
                     method="highs", options={"presolve": True,
                                             "dual_feasibility_tolerance": 1e-8,
                                             "primal_feasibility_tolerance": 1e-8})
    solve_seconds = time.perf_counter() - solve_start
    if not result.success:
        raise RuntimeError(f"Prototype optimization failed ({result.status}): {result.message}")
    solution = result.x
    d = model.data
    T, I = d.T, d.I
    effective = (model.capacity @ solution).reshape(T, I, 3)
    generation = (model.generation @ solution).reshape(T, I, 3)
    co2 = model.emissions @ solution
    annual_cost = (model.costs @ solution).reshape(T, len(COST_CATEGORIES))
    total_gen = generation.sum(axis=(1, 2))
    equality_residual = float(np.max(np.abs(model.A_eq @ solution - model.b_eq)))
    demand_violation = float(max(0, np.max(model.A_ub @ solution - model.b_ub)))
    summary: dict[str, Any] = {
        "case": "continuous_prototype", "status": "optimal",
        "solver": "scipy.optimize.linprog / HiGHS",
        "solver_message": result.message, "solver_status_code": int(result.status),
        "solver_iterations": int(result.nit),
        "runtime_versions": {"python": platform.python_version(), "numpy": np.__version__,
                             "scipy": scipy_version},
        "source_workbook": str(d.path),
        "source_sha256": hashlib.sha256(d.path.read_bytes()).hexdigest(),
        "horizon": T, "technologies": I, "variables": len(solution),
        "equality_constraints": model.A_eq.shape[0],
        "inequality_constraints": model.A_ub.shape[0],
        "objective": float(result.fun), "npv": float(model.npv @ solution),
        "latent_penalty_unweighted": float(model.latent @ solution),
        "latent_penalty": float(latent_factor * (model.latent @ solution)),
        "terminal_cost": float(model.terminal @ solution),
        "terminal_cost_weighted": float(1e-6 * (model.terminal @ solution)),
        "total_co2": float(co2.sum()), "total_generation": float(total_gen.sum()),
        "max_equality_residual": equality_residual,
        "max_demand_violation": demand_violation,
        "minimum_demand_margin_after_year0": float(np.min(total_gen[1:] - d.demand[1:T])),
        "build_seconds": model.build_seconds, "solve_seconds": solve_seconds,
        "elapsed_seconds": time.perf_counter() - started,
        "units": {"capacity": "GW", "generation": "TWh/year", "demand": "TWh/year",
                  "co2": "Mt CO2/year", "total_co2": "Mt CO2", "cost": "million USD (discounted)",
                  "objective": "billion USD (augmented discounted objective)"},
        "options": {"latent_factor": latent_factor, "no_retrofit": no_retrofit,
                    "no_delay_z": no_delay_z, "no_delay_x": no_delay_x},
        "preserved_source_quirks": [
            "Demand constraints apply to years 1..T-1, not year 0.",
            "Effective capacity is post-retirement/post-transition stock.",
            "Capital/fixed/variable modifiers use parent cost series even when base selectors are set.",
            "Retrofit latent penalty indexes service life by child kind.",
            "Retrofit terminal lost-sales coefficient omits capacity factor and uses default retirement tag."],
        "export_corrections": [
            "Existing age j=0 is included (old scripts incorrectly iterated 1..N).",
            "Generation sums Wgen+Zgen+Xgen; W is effective capacity, not generation.",
            "Workbook scaling already yields GW and TWh; demand/generation are not rescaled twice.",
            "Operating capacity uses W/Z/X, not pre-retirement w/z/x stocks."],
    }
    out = Path(output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    short = list(SHORT_NAMES) if I == len(SHORT_NAMES) else [f"tech_{i}" for i in range(I)]
    _write_csv(out / "capacity_by_tech.csv", ["year", *short],
               ([t, *effective[t].sum(axis=1)] for t in range(T)))
    _write_csv(out / "generation_by_tech.csv", ["year", *short],
               ([t, *generation[t].sum(axis=1)] for t in range(T)))
    _write_csv(out / "co2_emissions.csv", ["year", "co2_Mt"], enumerate(co2))
    _write_csv(out / "demand.csv", ["year", "demand_TWh", "generation_TWh", "margin_TWh", "constraint_enforced"],
               ([t, d.demand[t], total_gen[t], total_gen[t] - d.demand[t], t > 0] for t in range(T)))
    _write_csv(out / "costs.csv", ["year", *[f"{c}_million_USD" for c in COST_CATEGORIES], "total_million_USD"],
               ([t, *annual_cost[t], annual_cost[t].sum()] for t in range(T)))
    _write_csv(out / "assets_by_form.csv", ["year", "technology", "form", "capacity_GW", "generation_TWh"],
               ([t, short[i], form, effective[t, i, f], generation[t, i, f]]
                for t in range(T) for i in range(I)
                for f, form in enumerate(("existing", "retrofit", "new"))))
    _write_csv(out / "decisions.csv", ["variable", "year", "technology", "kind", "cohort", "value_GW"],
               ([group, key[0], short[key[1]], key[2] if len(key) == 4 or group == "xAlloc" else "",
                 key[-1] if len(key) == 4 or group == "uw" else "", solution[index]]
                for group in ("uw", "zTrans", "uz", "ux", "xAlloc")
                for key, index in model.variables[group].items() if solution[index] > 1e-8))
    _write_csv(out / "solution_variables.csv", ["variable", "index", "value"],
               ([group, str(key), float(solution[index])]
                for group, variables in model.variables.items()
                for key, index in variables.items()))
    _write_csv(out / "tech_names.csv", ["id", "technology", "name"],
               ([i, short[i], d.names[i]] for i in range(I)))
    (out / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    (out / "solver.log").write_text(
        f"Solver: {summary['solver']}\nStatus: {result.status}: {result.message}\n"
        f"Iterations: {result.nit}\nObjective: {result.fun:.15g}\n"
        f"Variables: {len(solution)}; equalities: {model.A_eq.shape[0]}; inequalities: {model.A_ub.shape[0]}\n"
        f"Maximum equality residual: {equality_residual:.12g}\n"
        f"Maximum demand violation: {demand_violation:.12g}\n"
        f"Build seconds: {model.build_seconds:.3f}; solve seconds: {solve_seconds:.3f}\n",
        encoding="utf-8")
    return summary
