"""Small sparse linear modeling layer backed by the native HiGHS Python API.

Free variables are the default, matching JuMP. Every row retains its original
constraint name so residuals can be traced back to the mathematical model.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import time
import os
from pathlib import Path
import shutil
import tempfile
import numpy as np


class Expr:
    __slots__ = ('terms', 'constant')
    __array_priority__ = 1000

    def __init__(self, terms=None, constant=0.0):
        self.terms = terms or {}
        self.constant = float(constant)

    def __add__(self, other):
        other = as_expr(other)
        terms = self.terms.copy()
        for col, coef in other.terms.items():
            terms[col] = terms.get(col, 0.0) + coef
        return Expr(terms, self.constant + other.constant)

    __radd__ = __add__

    def __neg__(self):
        return self * -1

    def __sub__(self, other):
        return self + (-as_expr(other))

    def __rsub__(self, other):
        return as_expr(other) + (-self)

    def __mul__(self, other):
        if isinstance(other, Expr):
            if other.terms and self.terms:
                raise TypeError('Nonlinear product in linear model')
            if other.terms:
                return other * self.constant
            other = other.constant
        factor = float(other)
        return Expr({j: c*factor for j, c in self.terms.items()}, self.constant*factor)

    __rmul__ = __mul__

    def __truediv__(self, other):
        return self * (1.0 / float(other))


class Var(Expr):
    __slots__ = ('index',)

    def __init__(self, index):
        super().__init__({index: 1.0})
        self.index = index


def as_expr(value):
    return value if isinstance(value, Expr) else Expr(constant=value)


def quicksum(items):
    terms = {}
    constant = 0.0
    for value in items:
        value = as_expr(value)
        constant += value.constant
        for col, coef in value.terms.items():
            terms[col] = terms.get(col, 0.0) + coef
    return Expr(terms, constant)


class IndexedVars(dict):
    def __getitem__(self, key):
        return super().__getitem__(key if isinstance(key, tuple) else (key,))

    def __setitem__(self, key, value):
        super().__setitem__(key if isinstance(key, tuple) else (key,), value)


@dataclass
class Row:
    expression: Expr
    lower: float
    upper: float
    name: str
    active: bool = True


class Model:
    def __init__(self):
        self.variables = {}
        self.column_names = []
        self.lower = []
        self.upper = []
        self.integer = []
        self.rows = []
        self.objective = Expr()
        self.solution = None
        self.solver = None

    def __getitem__(self, name):
        return self.variables[name]

    def var(self, name, key=(), lb=-math.inf, ub=math.inf, binary=False):
        if name not in self.variables:
            self.variables[name] = IndexedVars()
            setattr(self, name, self.variables[name])
        key = key if isinstance(key, tuple) else (key,)
        if key in self.variables[name]:
            raise ValueError(f'Duplicate variable {name}{key}')
        if binary:
            lb, ub = max(0.0, lb), min(1.0, ub)
        variable = Var(len(self.column_names))
        self.variables[name][key] = variable
        self.column_names.append(name + '[' + ','.join(map(str, key)) + ']')
        self.lower.append(float(lb))
        self.upper.append(float(ub))
        self.integer.append(bool(binary))
        return variable

    def add(self, lhs, sense, rhs=0.0, name=None):
        expr = as_expr(lhs) - rhs
        bound = -expr.constant
        expr.constant = 0.0
        if sense in ('==', '='):
            lower = upper = bound
        elif sense in ('<=', '≤'):
            lower, upper = -math.inf, bound
        elif sense in ('>=', '≥'):
            lower, upper = bound, math.inf
        else:
            raise ValueError(f'Unknown row sense: {sense}')
        self.rows.append(Row(expr, lower, upper, name or f'row_{len(self.rows)}'))
        return len(self.rows)-1

    def remove_group(self, name):
        for row in self.rows:
            if row.name == name or row.name.startswith(name + ':'):
                row.active = False

    def minimize(self, expr):
        self.objective = as_expr(expr)

    def value(self, expr):
        if self.solution is None:
            raise RuntimeError('No feasible solver solution available')
        expr = as_expr(expr)
        return float(expr.constant + sum(coef*self.solution[col] for col, coef in expr.terms.items()))

    def solve(self, time_limit=300.0, mip_gap=0.05, threads=4, log_file=None,
              model_file=None, starts=None, output=False):
        import highspy
        from scipy.sparse import csr_matrix
        active = [row for row in self.rows if row.active]
        rows, cols, data = [], [], []
        for i, row in enumerate(active):
            for j, coef in row.expression.terms.items():
                if coef:
                    rows.append(i); cols.append(j); data.append(coef)
        matrix = csr_matrix((data, (rows, cols)), shape=(len(active), len(self.lower)))
        model = highspy.HighsLp()
        model.num_col_ = len(self.lower)
        model.num_row_ = len(active)
        cost = np.zeros(model.num_col_)
        for col, coef in self.objective.terms.items():
            cost[col] = coef
        model.col_cost_ = cost
        model.col_lower_ = self.lower
        model.col_upper_ = self.upper
        model.row_lower_ = [r.lower for r in active]
        model.row_upper_ = [r.upper for r in active]
        model.offset_ = self.objective.constant
        model.col_names_ = self.column_names
        model.row_names_ = [r.name for r in active]
        model.a_matrix_.format_ = highspy.MatrixFormat.kRowwise
        model.a_matrix_.start_ = matrix.indptr
        model.a_matrix_.index_ = matrix.indices
        model.a_matrix_.value_ = matrix.data
        model.integrality_ = [highspy.HighsVarType.kInteger if integer else
                              highspy.HighsVarType.kContinuous for integer in self.integer]
        h = highspy.Highs()
        # HiGHS' Windows C++ file API cannot open every Unicode path. Stage
        # native solver files under an ASCII temporary directory; Python copies
        # them to the caller's chosen (possibly Chinese) directory.
        file_stage = tempfile.TemporaryDirectory(prefix='stream_highs_',dir=Path.cwd())
        stage = Path(file_stage.name)
        native_stage = Path(os.path.relpath(stage))
        # A fresh scheduler is necessary when another case used different threads.
        highspy.Highs.resetGlobalScheduler(True)
        for key, value in dict(time_limit=float(time_limit), mip_rel_gap=float(mip_gap),
                               threads=int(threads), output_flag=bool(output), random_seed=0).items():
            h.setOptionValue(key, value)
        if log_file:
            h.setOptionValue('log_file', str(native_stage/'solver.log'))
            h.setOptionValue('output_flag', True)
            h.setOptionValue('log_to_console', bool(output))
        status = h.passModel(model)
        if status == highspy.HighsStatus.kError:
            raise RuntimeError('HiGHS rejected the generated model')
        if model_file:
            if h.writeModel(str(native_stage/'model.mps'))==highspy.HighsStatus.kError:
                raise RuntimeError(f'Cannot write solver model file: {model_file}')
            shutil.copy2(stage/'model.mps',model_file)
        if starts:
            indices = np.asarray(list(starts), dtype=np.int32)
            values = np.asarray(list(starts.values()), dtype=float)
            h.setSolution(len(indices), indices, values)
        start = time.perf_counter()
        h.run()
        seconds = time.perf_counter()-start
        if log_file:
            h.setOptionValue('log_file','')
            shutil.copy2(stage/'solver.log',log_file)
        file_stage.cleanup()
        info = h.getInfo()
        status_name = h.modelStatusToString(h.getModelStatus())
        sol = h.getSolution()
        if not sol.value_valid or info.primal_solution_status != highspy.SolutionStatus.kSolutionStatusFeasible:
            raise RuntimeError(f'HiGHS {status_name}; no feasible incumbent, nothing exported')
        self.solution = np.asarray(sol.col_value)
        self.solver = h
        activity = matrix @ self.solution
        lower = np.asarray(model.row_lower_)
        upper = np.asarray(model.row_upper_)
        violations = np.maximum(np.maximum(lower - activity, activity - upper), 0.0)
        scale = np.maximum(1.0, np.maximum(np.abs(np.where(np.isfinite(lower),lower,0)),
                                         np.abs(np.where(np.isfinite(upper),upper,0))))
        bound_v = max(float(np.max(np.maximum(np.asarray(self.lower)-self.solution, 0), initial=0)),
                      float(np.max(np.maximum(self.solution-np.asarray(self.upper), 0), initial=0)))
        integer_x = self.solution[np.asarray(self.integer)]
        integrality = float(np.max(np.abs(integer_x-np.round(integer_x)), initial=0))
        worst = np.argsort(violations)[-10:][::-1]
        return dict(status=status_name, objective=self.value(self.objective),
                    primal_bound=info.objective_function_value,
                    dual_bound=info.mip_dual_bound if any(self.integer) else info.objective_function_value,
                    mip_gap=info.mip_gap if any(self.integer) else 0.0,
                    solve_seconds=seconds, variables=len(self.lower), constraints=len(active),
                    binary_variables=sum(self.integer), nonzeros=matrix.nnz,
                    max_constraint_violation=float(np.max(violations, initial=0)),
                    max_relative_constraint_violation=float(np.max(violations/scale, initial=0)),
                    max_bound_violation=bound_v, max_integrality_violation=integrality,
                    worst_rows=[dict(name=active[i].name, violation=float(violations[i])) for i in worst],
                    highs_version=h.version())
