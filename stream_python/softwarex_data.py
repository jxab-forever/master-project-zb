"""Excel parameters for the SoftwareX STREAM model.

Port of Argonne STREAM's ``kern.jl:read_params`` (BSD-3-Clause).
Excel flattening follows the reader, including its overlapping onsite columns.
No Julia installation or source checkout is needed at runtime.
"""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import numbers
import numpy as np
from openpyxl import load_workbook


class JArray:
    """Numpy storage with Julia's 1-based, column-major indexing.

    Slice stops are inclusive, as in Julia ranges. ``a[:]`` is a linear
    column-major vector; ``a[:, 2]`` selects a column. Use ``a.data`` or
    ``np.asarray(a)`` for ordinary zero-based numpy operations.
    """
    def __init__(self, data, dtype=None):
        self.data = np.asarray(data, dtype=dtype)

    @property
    def shape(self):
        return self.data.shape

    @property
    def ndim(self):
        return self.data.ndim

    @property
    def size(self):
        return self.data.size

    @property
    def dtype(self):
        return self.data.dtype

    def __array__(self, dtype=None, copy=None):
        a = np.asarray(self.data, dtype=dtype)
        return a.copy() if copy else a

    def __len__(self):
        return self.size

    def __iter__(self):
        return iter(self.data.ravel(order="F"))

    @staticmethod
    def _index(key, size):
        if isinstance(key, numbers.Integral):
            if not 1 <= key <= size:
                raise IndexError(f"Julia index {key} outside 1:{size}")
            return int(key) - 1
        if isinstance(key, slice):
            step = 1 if key.step is None else key.step
            start = (1 if step > 0 else size) if key.start is None else key.start
            stop = (size if step > 0 else 1) if key.stop is None else key.stop
            # Explicit vectors also handle a backwards range ending at 1.
            return np.array([JArray._index(k, size) for k in range(start, stop + (1 if step > 0 else -1), step)], dtype=int)
        if isinstance(key, (range, list, np.ndarray)):
            return np.array([JArray._index(k, size) for k in key], dtype=int)
        raise TypeError(f"Unsupported Julia index: {key!r}")

    def _selection(self, key):
        if not isinstance(key, tuple):
            return self.data.reshape(-1, order="F"), self._index(key, self.size)
        a = self.data
        if len(key) > a.ndim:
            a = a.reshape(a.shape + (1,) * (len(key) - a.ndim))
        if len(key) < a.ndim:
            # Julia collapses omitted trailing dimensions into the last index.
            a = a.reshape(a.shape[:len(key)-1] + (int(np.prod(a.shape[len(key)-1:])),), order="F")
        idx = [self._index(k, dim) for k, dim in zip(key, a.shape)]
        # Cartesian selection (numpy's advanced indexing is pairwise by default).
        vectors = [i for i in idx if isinstance(i, np.ndarray)]
        if len(vectors) > 1:
            meshes = iter(np.ix_(*vectors))
            idx = [next(meshes) if isinstance(i, np.ndarray) else i for i in idx]
        return a, tuple(idx)

    def __getitem__(self, key):
        a, idx = self._selection(key)
        value = a[idx]
        return JArray(value) if isinstance(value, np.ndarray) else value.item() if isinstance(value, np.generic) else value

    def __setitem__(self, key, value):
        # A flat view of non-Fortran contiguous data may be a copy.
        if not isinstance(key, tuple):
            idx = self._index(key, self.size)
            coords = np.unravel_index(idx, self.shape, order="F")
            self.data[coords] = np.asarray(value)
        else:
            a, idx = self._selection(key)
            a[idx] = np.asarray(value)

    def copy(self):
        return JArray(self.data.copy())

    def __repr__(self):
        return f"JArray({self.data!r})"


class JList(list):
    """A 1-based sequence for ragged arrays and link tuples."""
    def __getitem__(self, key):
        if isinstance(key, slice):
            return JList(list.__getitem__(self, int(k)) for k in JArray._index(key, len(self)))
        return list.__getitem__(self, JArray._index(key, len(self)))

    def __setitem__(self, key, value):
        list.__setitem__(self, JArray._index(key, len(self)), value)


# (parameter attribute, worksheet, shape dimensions, scalar type).
# Generated once from the direct rectangular reads in kern.jl.
_DIRECT = [('n_rfu', 'n_rfu', ('n_loc',), 'int'),
 ('n_nfu', 'n_nfu', ('n_loc',), 'int'),
 ('c0', 'c0', ('n_loc',), 'float'),
 ('e_C', 'e_C', ('n_loc',), 'float'),
 ('e_c_ub', 'e_c_ub', ('n_loc',), 'float'),
 ('e_loanFact', 'e_loanFact', ('n_loc',), 'float'),
 ('e_l_ub', 'e_l_ub', ('n_loc',), 'float'),
 ('e_Ann', 'e_Ann', ('n_loc',), 'float'),
 ('e_ann_ub', 'e_ann_ub', ('n_loc',), 'float'),
 ('e_ladd_ub', 'e_ladd_ub', ('n_loc',), 'float'),
 ('e_loan_ub', 'e_loan_ub', ('n_loc',), 'float'),
 ('e_pay_ub', 'e_pay_ub', ('n_loc',), 'float'),
 ('r_filter', 'r_filter', ('n_loc', 'n_rtft'), 'bool'),
 ('r_c_C', 'r_c_C', ('n_loc', 'n_rtft'), 'float'),
 ('r_rhs_C', 'r_rhs_C', ('n_loc', 'n_rtft'), 'float'),
 ('r_cp_ub', 'r_cp_ub', ('n_loc', 'n_node'), 'float'),
 ('r_cpb_ub', 'r_cpb_ub', ('n_loc',), 'float'),
 ('r_eh_ub', 'r_eh_ub', ('n_loc', 'n_node'), 'float'),
 ('r_c_Hfac', 'r_c_Hfac', ('n_loc', 'n_rtft'), 'float'),
 ('r_c_Helec', 'r_c_Helec', ('n_loc', 'n_rtft'), 'float'),
 ('r_c_Ufac', 'r_c_Ufac', ('n_loc', 'n_rtft'), 'float'),
 ('r_u_ub', 'r_u_ub', ('n_loc', 'n_node'), 'float'),
 ('r_cpe_ub', 'r_cpe_ub', ('n_loc', 'n_node'), 'float'),
 ('r_fu_e_ub', 'r_fu_e_ub', ('n_loc', 'n_node'), 'float'),
 ('r_u_fu_e_ub', 'r_u_fu_e_ub', ('n_loc', 'n_node'), 'float'),
 ('r_ep0_ub', 'r_ep0_ub', ('n_loc',), 'float'),
 ('r_chi', 'r_chi', ('n_loc', 'n_rtft'), 'float'),
 ('r_ep1ge_ub', 'r_ep1ge_ub', ('n_loc',), 'float'),
 ('r_sigma', 'r_sigma', ('n_loc', 'n_rtft'), 'float'),
 ('r_ep1gce_ub', 'r_ep1gce_ub', ('n_loc',), 'float'),
 ('r_ep1gcs_ub', 'r_ep1gcs_ub', ('n_loc',), 'float'),
 ('r_c_fOnm', 'r_c_fOnm', ('n_loc', 'n_rtft'), 'float'),
 ('r_rhs_fOnm', 'r_rhs_fOnm', ('n_loc', 'n_rtft'), 'float'),
 ('r_cfonm_ub', 'r_cfonm_ub', ('n_loc', '1'), 'float'),
 ('r_c_vOnm', 'r_c_vOnm', ('n_loc', 'n_rtft'), 'float'),
 ('r_rhs_vOnm', 'r_rhs_vOnm', ('n_loc', 'n_rtft'), 'float'),
 ('r_cvonm_ub', 'r_cvonm_ub', ('n_loc', '1'), 'float'),
 ('r_e_c_ub', 'r_e_c_ub', ('n_loc',), 'float'),
 ('r_l0_ub', 'r_l0_ub', ('n_loc', 'n_rtft'), 'float'),
 ('r_le_ub', 'r_le_ub', ('n_loc', 'n_rtft'), 'float'),
 ('r_ann0_bM', 'r_ann0_bM', ('n_loc', '1'), 'float'),
 ('r_anne_bM', 'r_anne_bM', ('n_loc', '1'), 'float'),
 ('r_l0add_bM', 'r_l0add_bM', ('n_loc', '1'), 'float'),
 ('r_leadd_bM', 'r_leadd_bM', ('n_loc', '1'), 'float'),
 ('r_loan_ub', 'r_loan_ub', ('n_loc', '1'), 'float'),
 ('r_pay0_bM', 'r_pay0_ub', ('n_loc', '1'), 'float'),
 ('r_paye_bM', 'r_paye_ub', ('n_loc', '1'), 'float'),
 ('r_c_upsein_rate', 'r_c_upsein_rate', ('n_loc', 'n_mat'), 'float'),
 ('r_ups_e_mt_in_ub', 'r_ups_e_mt_in_ub', ('n_loc',), 'float'),
 ('n_filter', 'n_filter', ('n_loc', 'n_new'), 'bool'),
 ('n_cp_bM', 'n_cp_bM', ('n_loc', 'n_node'), 'float'),
 ('n_c0_bM', 'n_c0_bM', ('n_loc',), 'float'),
 ('n_c0_lo', 'n_c0_lo', ('n_loc', 'n_new'), 'float'),
 ('n_loanFact', 'n_loanFact', ('n_loc', 'n_new'), 'float'),
 ('n_l_bM', 'n_l_bM', ('n_loc', '1'), 'float'),
 ('n_Ann', 'n_Ann', ('n_loc', 'n_new'), 'float'),
 ('n_ann_bM', 'n_ann_bM', ('n_loc', '1'), 'float'),
 ('n_ladd_bM', 'n_ladd_bM', ('n_loc', '1'), 'float'),
 ('n_loan_bM', 'n_loan_bM', ('n_loc', '1'), 'float'),
 ('n_pay_bM', 'n_pay_bM', ('n_loc', '1'), 'float'),
 ('n_eh_ub', 'n_eh_ub', ('n_loc', 'n_node'), 'float'),
 ('n_c_Hfac', 'n_c_Hfac', ('n_loc', 'n_new'), 'float'),
 ('n_c_Helec', 'n_c_Helec', ('n_loc', 'n_new'), 'float'),
 ('n_u_ub', 'n_u_ub', ('n_loc', 'n_node'), 'float'),
 ('n_cpe_ub', 'n_cpe_ub', ('n_loc', 'n_node'), 'float'),
 ('n_ep0_bM', 'n_ep0_bM', ('n_loc',), 'float'),
 ('n_chi', 'n_chi', ('n_loc', 'n_new'), 'float'),
 ('n_ep1ge_bM', 'n_ep1ge_bM', ('n_loc',), 'float'),
 ('n_sigma', 'n_sigma', ('n_loc', 'n_new'), 'float'),
 ('n_ep1gce_bM', 'n_ep1gce_bM', ('n_loc',), 'float'),
 ('n_ep1gcs_bM', 'n_ep1gcs_bM', ('n_loc',), 'float'),
 ('n_c_fOnm', 'n_c_fOnm', ('n_loc', 'n_new'), 'float'),
 ('n_rhs_fOnm', 'n_rhs_fOnm', ('n_loc', 'n_new'), 'float'),
 ('n_cfonm_bM', 'n_cfonm_bM', ('n_loc',), 'float'),
 ('n_c_vOnm', 'n_c_vOnm', ('n_loc', 'n_new'), 'float'),
 ('n_rhs_vOnm', 'n_rhs_vOnm', ('n_loc', 'n_new'), 'float'),
 ('n_cvonm_bM', 'n_cvonm_bM', ('n_loc',), 'float'),
 ('n_c_upsein_rate', 'n_c_upsein_rate', ('n_loc', 'n_mat'), 'float'),
 ('n_ups_e_mt_in_ub', 'n_ups_e_mt_in_ub', ('n_loc',), 'float'),
 ('c_cts_cost', 'c_cts_cost', ('n_loc',), 'float'),
 ('c_xin_cost', 'c_xin_cost', ('n_loc', 'n_mat'), 'float'),
 ('o_cp_ub', 'o_cp_ub', ('n_loc', 'n_node'), 'float'),
 ('o_u_ub', 'o_u_ub', ('n_loc', 'n_node'), 'float'),
 ('o_ep0_bM', 'o_ep0_bM', ('n_loc',), 'float'),
 ('o_ep1ge_bM', 'o_ep1ge_bM', ('n_loc',), 'float'),
 ('o_ep1gce_bM', 'o_ep1gce_bM', ('n_loc',), 'float'),
 ('o_ep1gcs_bM', 'o_ep1gcs_bM', ('n_loc',), 'float'),
 ('o_ups_e_mt_in_ub', 'o_ups_e_mt_in_ub', ('n_loc',), 'float'),
 ('o_pay_bM', 'o_pay_bM', ('n_loc', '1'), 'float'),
 ('o_cfonm_bM', 'o_cfonm_bM', ('n_loc', '1'), 'float'),
 ('o_cvonm_bM', 'o_cvonm_bM', ('n_loc', '1'), 'float'),
 ('o_fstck_ub', 'o_fstck_ub', ('n_loc', 'n_fstck'), 'float'),
 ('t_ret_c_bM', 't_ret_c_bM', ('n_loc', '1'), 'float'),
 ('t_loan_bM', 't_loan_bM', ('n_loc', '1'), 'float'),
 ('r_loan0', 'r_loan0', ('n_loc',), 'float'),
 ('node_mat', 'node_mat', ('n_node', 'n_mat'), 'bool'),
 ('skip_mb', 'skip_mb', ('n_node', 'n_mat'), 'bool'),
 ('ckey', 'ckey', ('n_node',), 'int'),
 ('nd_en_fltr', 'nd_en_fltr', ('n_node',), 'bool'),
 ('nd_em_fltr', 'nd_em_fltr', ('n_node',), 'bool'),
 ('r_Kkey_j', 'r_Kkey_j', ('n_node', 'n_rtft'), 'float'),
 ('n_Kkey_j', 'n_Kkey_j', ('n_node', 'n_new'), 'float'),
 ('min_cpr', 'min_cpr', ('n_loc',), 'float')]


def read_params(path: str | Path) -> SimpleNamespace:
    """Read all original ``params`` fields with the Julia reader's layout."""
    wb = load_workbook(path, data_only=True, read_only=False)
    try:
        first = wb.worksheets[0]
        scalar_names = "n_periods n_subperiods n_location n_rtft n_new n_fstck n_node n_mat n_link yr_subperiod y0 x_ub interest sf_cap sf_cash sf_heat sf_elec sf_em key_node initial_co2".split()
        values = {}
        for row, name in enumerate(scalar_names, 2):
            value = first.cell(row, 2).value
            values[name] = int(value) if row <= 12 or name == "key_node" else float(value)
        p = SimpleNamespace(**values)
        dims = {**values, "n_loc": p.n_location, "1": 1}

        def cell(sheet, row, col):
            v = wb[sheet].cell(row, col).value
            if v is None:
                raise ValueError(f"Missing parameter value at {sheet}!R{row}C{col}")
            return v

        def rectangle(sheet, rows, cols, dtype=float):
            return np.array([[cell(sheet, r+2, c+1) for c in range(cols)] for r in range(rows)], dtype=dtype)

        def mapped(name, shape, coordinates, sheet=None):
            data = np.empty(shape, dtype=float)
            for idx in np.ndindex(shape):
                row, col = coordinates(*idx)
                data[idx] = cell(sheet or name, row+2, col+1)
            setattr(p, name, JArray(data))

        for attr, sheet, shape, dtype in _DIRECT:
            shape = tuple(dims[d] for d in shape)
            a = rectangle(sheet, shape[0], shape[1] if len(shape) == 2 else 1, {"bool": bool, "int": int, "float": float}[dtype])
            setattr(p, attr, JArray(a.reshape(shape)))
        I, J, L, N, M, F = p.n_periods, p.n_subperiods, p.n_location, p.n_node, p.n_mat, p.n_fstck
        for prefix, K, fuel in (("r", p.n_rtft, max(p.n_rfu)), ("n", p.n_new, max(p.n_nfu))):
            for suffix in ("c_H", "rhs_H", "c_U", "rhs_U", "c_cpe", "rhs_cpe"):
                mapped(f"{prefix}_{suffix}", (L,K,N), lambda l,k,n: (l,k+K*n))
            for suffix in ("c_F", "rhs_F", "c_Fe", "c_Fgenf", "c_Hr"):
                mapped(f"{prefix}_{suffix}", (fuel,L,K,N), lambda f,l,k,n: (f,l+L*k+L*K*n))
            for suffix in ("ehf_ub", "u_ehf_ub"):
                mapped(f"{prefix}_{suffix}", (fuel,L,N), lambda f,l,n: (f,l+L*n))
            mapped(f"{prefix}_c_UonSite", (I,J,L,K,N), lambda i,j,l,k,n: (j+J*i,k+K*l+L*n))
            for suffix in ("c_Fstck", "rhs_Fstck", "fstck_ub"):
                mapped(f"{prefix}_{suffix}", (L,K,F), lambda l,k,f: (l,f+F*k))
            mapped(f"{prefix}_Kmb", (N,M,M,K), lambda n,m,c,k: (m,k+K*c+K*M*n))
            for suffix in ("x_in_ub", "x_out_ub"):
                mapped(f"{prefix}_{suffix}", (L,N,M), lambda l,n,c: (l,c+M*n))
        # These three sheets are loaded via loops instead of rectangular reads.
        p.n_c_Ufac = JArray(rectangle("n_c_Ufac", L, p.n_new))
        for name in ("n_fu_e_ub", "n_u_fu_e_ub"):
            setattr(p, name, JArray(rectangle(name, L, N)))
        for name in ("r_loanFact", "r_Ann"):
            mapped(name, (I,J,L,p.n_rtft), lambda i,j,l,k: (j+J*i,k+p.n_rtft*l))
        for name in ("c_u_cost", "GcI"):
            mapped(name, (I,J,L), lambda i,j,l: (j+J*i,l))
        for name, fuel in (("c_r_ehf_cost", max(p.n_rfu)), ("c_n_ehf_cost", max(p.n_nfu))):
            mapped(name, (I,J,L,fuel), lambda i,j,l,f: (j+J*i,f+fuel*l))
        mapped("o_ehf_ub", (max(p.n_rfu),L,N), lambda f,l,n: (f,l+L*n))
        for name in ("o_x_in_ub", "o_x_out_ub"):
            mapped(name, (L,N,M), lambda l,n,c: (l,c+M*n))
        for name in ("discount", "demand", "co2_budget", "elec_budget"):
            mapped(name, (I,J), lambda i,j: (j+J*i,0))
        p.o_cpe_bM = float(cell("o_cpe_bM", 2, 1))
        for name in ("input_mat", "output_mat"):
            setattr(p, name, JList(JList(int(cell(name,r+3,n+1)) for r in range(int(cell(name,2,n+1)))) for n in range(N)))
        p.links_list = JList(JList(int(cell("links_list",r+2,c+1)) for c in range(3)) for r in range(p.n_link))
        return p
    finally:
        wb.close()


def sets(p: SimpleNamespace) -> SimpleNamespace:
    """Return the original 1-based model index sets."""
    names = {"P":"n_periods", "P2":"n_subperiods", "L":"n_location", "Kr":"n_rtft", "Kn":"n_new", "Nf":"n_fstck", "Nd":"n_node", "Mt":"n_mat", "Ln":"n_link"}
    result = {name:range(1, getattr(p, dim)+1) for name, dim in names.items()}
    result["Fu_r"] = JList(range(1,int(n)+1) for n in p.n_rfu)
    result["Fu_n"] = JList(range(1,int(n)+1) for n in p.n_nfu)
    return SimpleNamespace(**result)
