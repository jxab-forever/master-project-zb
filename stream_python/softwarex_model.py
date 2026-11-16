"""SoftwareX STRE3AM model, translated equation by equation from bm.jl.

Copyright (C) 2023, UChicago Argonne, LLC. All Rights Reserved.
Software: STRE3AM, by Argonne National Laboratory. BSD-3-Clause.
See LICENSE for the complete original license and disclaimer.
Generated with tools/port_softwarex.py; original line numbers accompany equations.
All indices retain the original 1-based convention. No LP snapshots or Julia calls.
"""
from math import inf
import numpy as np
from .linear import Model, quicksum

sum = quicksum

class IndexedVars(dict):
    def __getitem__(self, key):
        return super().__getitem__(key if isinstance(key, tuple) else (key,))

def _in(a, b): return a in b
def _ne(a, b): return a != b
def keys(a): return a.keys()
def _range(a, b): return range(int(a), int(b) + 1)
def _vec(a): return {i + 1: x for i, x in enumerate(np.asarray(a).ravel(order="F"))}
def _maximum(a, dims=None):
    a = np.asarray(a)
    if dims is None: return float(np.max(a))
    return {i + 1: float(x) for i, x in enumerate(np.max(a, axis=dims-1).ravel(order="F"))}

def gen_lp_window(p, s):
    slots = [(i, j) for i in s.P for j in s.P2]
    loan = {slot: slots[k:] for k, slot in enumerate(slots)}
    pay = {slot: slots[k:k + 30 // p.yr_subperiod] for k, slot in enumerate(slots)}
    return loan, pay

def f_xin(m, i, j, l, k, n, c, input_mat, mode="r"):
    return getattr(m, mode + "_x_in_d_")[i,j,l,k,n,c] if c in input_mat[n] else 0.0

def f_xout(m, i, j, l, k, n, c, output_mat, mode="n"):
    return getattr(m, mode + "_x_out_d_")[i,j,l,k,n,c] if c in output_mat[n] else 0.0


def createBlockMod(index_p, index_l, p, s):
    P = s.P
    P2 = s.P2
    L = s.L
    Kr = s.Kr
    Kn = s.Kn
    Fu_r = s.Fu_r
    Fu_n = s.Fu_n
    Nf = s.Nf
    Nd = s.Nd
    Mt = s.Mt
    Ln = s.Ln
    n_periods = p.n_periods
    n_subperiods = p.n_subperiods
    yr_subperiod = p.yr_subperiod
    loan_subperiods = 30 // yr_subperiod
    loan_periods = loan_subperiods // n_subperiods
    remain_loan_years = 30 - loan_subperiods * yr_subperiod
    remain_subp_loan = loan_subperiods - loan_periods * n_subperiods
    remain_sub_loan = remain_subp_loan + int(remain_loan_years > 0)
    sT = 1
    sF = 2
    fP = min(P)
    fP2 = min(P2)
    fKr = min(Kr)
    fKn = min(Kn)
    r_filter = p.r_filter
    n_filter = p.n_filter
    nd_en_fltr = p.nd_en_fltr
    nd_em_fltr = p.nd_em_fltr
    node_mat = p.node_mat
    skip_mb = p.skip_mb
    links_list = p.links_list
    input_mat = p.input_mat
    output_mat = p.output_mat
    ckey = p.ckey
    key_node = p.key_node
    loan_window, pay_window = gen_lp_window(p, s)
    m = Model()
    P = [index_p] if isinstance(index_p, int) else index_p
    L = [index_l] if isinstance(index_l, int) else index_l
    # bm.jl:217 — variable y_o
    y_o = m.y_o = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                y_o[(i, j, l)] = m.var('y_o', (i, j, l), lb=-inf, ub=1.0, binary=True)
    # bm.jl:220 — variable y_r
    y_r = m.y_r = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        y_r[(i, j, l, k)] = m.var('y_r', (i, j, l, k), lb=-inf, ub=1.0, binary=True)
    # bm.jl:228 — variable y_e
    y_e = m.y_e = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                y_e[(i, j, l)] = m.var('y_e', (i, j, l), lb=-inf, ub=1.0, binary=True)
    # bm.jl:230 — variable e_c_d_
    e_c_d_ = m.e_c_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for _state in (sT, sF):
                    e_c_d_[(i, j, l, _state)] = m.var('e_c_d_', (i, j, l, _state), lb=0.0, ub=inf, binary=False)
    # bm.jl:232 — variable e_c
    e_c = m.e_c = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                e_c[(i, j, l)] = m.var('e_c', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:234 — variable e_l
    e_l = m.e_l = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                e_l[(i, j, l)] = m.var('e_l', (i, j, l), lb=0.0, ub=p.e_l_ub[l], binary=False)
    # bm.jl:236 — variable x
    x = m.x = IndexedVars()
    for i in P:
        for l in L:
            x[(i, l)] = m.var('x', (i, l), lb=0.0, ub=p.x_ub, binary=False)
    # bm.jl:238 — variable x_d_
    x_d_ = m.x_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for _state in (sT, sF):
                    x_d_[(i, j, l, _state)] = m.var('x_d_', (i, j, l, _state), lb=0.0, ub=inf, binary=False)
    # bm.jl:243 — variable r_cp
    r_cp = m.r_cp = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    r_cp[(i, j, l, n)] = m.var('r_cp', (i, j, l, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:245 — variable r_cp_d_
    r_cp_d_ = m.r_cp_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]:
                            r_cp_d_[(i, j, l, k, n)] = m.var('r_cp_d_', (i, j, l, k, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:246 — variable cpb
    cpb = m.cpb = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                cpb[(i, j, l)] = m.var('cpb', (i, j, l), lb=0.0, ub=inf, binary=False)
    # bm.jl:248 — variable r_cpb_d_
    r_cpb_d_ = m.r_cpb_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        r_cpb_d_[(i, j, l, k)] = m.var('r_cpb_d_', (i, j, l, k), lb=0.0, ub=inf, binary=False)
    # bm.jl:251 — variable r_eh
    r_eh = m.r_eh = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        r_eh[(i, j, l, n)] = m.var('r_eh', (i, j, l, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:253 — variable r_eh_d_
    r_eh_d_ = m.r_eh_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]  and  nd_en_fltr[n]:
                            r_eh_d_[(i, j, l, k, n)] = m.var('r_eh_d_', (i, j, l, k, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:256 — variable r_ehf
    r_ehf = m.r_ehf = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for f in Fu_r[l]:
                    for n in Nd:
                        if nd_en_fltr[n]:
                            r_ehf[(i, j, l, f, n)] = m.var('r_ehf', (i, j, l, f, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:258 — variable r_ehf_d_
    r_ehf_d_ = m.r_ehf_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for f in Fu_r[l]:
                        for n in Nd:
                            if r_filter[l, k]  and  nd_en_fltr[n]:
                                r_ehf_d_[(i, j, l, k, f, n)] = m.var('r_ehf_d_', (i, j, l, k, f, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:261 — variable r_u
    r_u = m.r_u = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        r_u[(i, j, l, n)] = m.var('r_u', (i, j, l, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:263 — variable r_u_d_
    r_u_d_ = m.r_u_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]  and  nd_en_fltr[n]:
                            r_u_d_[(i, j, l, k, n)] = m.var('r_u_d_', (i, j, l, k, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:266 — variable r_u_onsite
    r_u_onsite = m.r_u_onsite = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        r_u_onsite[(i, j, l, n)] = m.var('r_u_onsite', (i, j, l, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:268 — variable r_u_onsite_d_
    r_u_onsite_d_ = m.r_u_onsite_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]  and  nd_en_fltr[n]:
                            r_u_onsite_d_[(i, j, l, k, n)] = m.var('r_u_onsite_d_', (i, j, l, k, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:271 — variable r_u_ehf
    r_u_ehf = m.r_u_ehf = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for f in Fu_r[l]:
                    for n in Nd:
                        if nd_en_fltr[n]:
                            r_u_ehf[(i, j, l, f, n)] = m.var('r_u_ehf', (i, j, l, f, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:273 — variable r_u_ehf_d_
    r_u_ehf_d_ = m.r_u_ehf_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for f in Fu_r[l]:
                        for n in Nd:
                            if r_filter[l, k]  and  nd_en_fltr[n]:
                                r_u_ehf_d_[(i, j, l, k, f, n)] = m.var('r_u_ehf_d_', (i, j, l, k, f, n), lb=0.0, ub=inf, binary=False)
    # bm.jl:276 — variable r_cpe
    r_cpe = m.r_cpe = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_em_fltr[n]:
                        r_cpe[(i, j, l, n)] = m.var('r_cpe', (i, j, l, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:278 — variable r_cpe_d_
    r_cpe_d_ = m.r_cpe_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]  and  nd_em_fltr[n]:
                            r_cpe_d_[(i, j, l, k, n)] = m.var('r_cpe_d_', (i, j, l, k, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:281 — variable r_fu_e
    r_fu_e = m.r_fu_e = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        r_fu_e[(i, j, l, n)] = m.var('r_fu_e', (i, j, l, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:283 — variable r_fu_e_d_
    r_fu_e_d_ = m.r_fu_e_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]  and  nd_en_fltr[n]:
                            r_fu_e_d_[(i, j, l, k, n)] = m.var('r_fu_e_d_', (i, j, l, k, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:286 — variable r_u_fu_e
    r_u_fu_e = m.r_u_fu_e = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        r_u_fu_e[(i, j, l, n)] = m.var('r_u_fu_e', (i, j, l, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:288 — variable r_u_fu_e_d_
    r_u_fu_e_d_ = m.r_u_fu_e_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]  and  nd_en_fltr[n]:
                            r_u_fu_e_d_[(i, j, l, k, n)] = m.var('r_u_fu_e_d_', (i, j, l, k, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:291 — variable r_ep0
    r_ep0 = m.r_ep0 = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                r_ep0[(i, j, l)] = m.var('r_ep0', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:293 — variable r_ep0_d_
    r_ep0_d_ = m.r_ep0_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        r_ep0_d_[(i, j, l, k)] = m.var('r_ep0_d_', (i, j, l, k), lb=-inf, ub=inf, binary=False)
    # bm.jl:295 — variable r_ep1ge
    r_ep1ge = m.r_ep1ge = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                r_ep1ge[(i, j, l)] = m.var('r_ep1ge', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:297 — variable r_ep1ge_d_
    r_ep1ge_d_ = m.r_ep1ge_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        r_ep1ge_d_[(i, j, l, k)] = m.var('r_ep1ge_d_', (i, j, l, k), lb=-inf, ub=inf, binary=False)
    # bm.jl:299 — variable r_ep1gce
    r_ep1gce = m.r_ep1gce = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                r_ep1gce[(i, j, l)] = m.var('r_ep1gce', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:301 — variable r_ep1gce_d_
    r_ep1gce_d_ = m.r_ep1gce_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        r_ep1gce_d_[(i, j, l, k)] = m.var('r_ep1gce_d_', (i, j, l, k), lb=-inf, ub=inf, binary=False)
    # bm.jl:303 — variable r_ep1gcs
    r_ep1gcs = m.r_ep1gcs = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                r_ep1gcs[(i, j, l)] = m.var('r_ep1gcs', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:305 — variable r_ep1gcs_d_
    r_ep1gcs_d_ = m.r_ep1gcs_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        r_ep1gcs_d_[(i, j, l, k)] = m.var('r_ep1gcs_d_', (i, j, l, k), lb=-inf, ub=inf, binary=False)
    # bm.jl:307 — variable r_fstck
    r_fstck = m.r_fstck = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for f in Nf:
                    r_fstck[(i, j, l, f)] = m.var('r_fstck', (i, j, l, f), lb=-inf, ub=inf, binary=False)
    # bm.jl:309 — variable r_fstck_d_
    r_fstck_d_ = m.r_fstck_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for f in Nf:
                        if r_filter[l, k]:
                            r_fstck_d_[(i, j, l, k, f)] = m.var('r_fstck_d_', (i, j, l, k, f), lb=0.0, ub=inf, binary=False)
    # bm.jl:311 — variable r_cfonm_d_
    r_cfonm_d_ = m.r_cfonm_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        r_cfonm_d_[(i, j, l, k)] = m.var('r_cfonm_d_', (i, j, l, k), lb=-inf, ub=inf, binary=False)
    # bm.jl:313 — variable r_cfonm
    r_cfonm = m.r_cfonm = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                r_cfonm[(i, j, l)] = m.var('r_cfonm', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:316 — variable r_cvonm_d_
    r_cvonm_d_ = m.r_cvonm_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        r_cvonm_d_[(i, j, l, k)] = m.var('r_cvonm_d_', (i, j, l, k), lb=-inf, ub=inf, binary=False)
    # bm.jl:318 — variable r_cvonm
    r_cvonm = m.r_cvonm = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                r_cvonm[(i, j, l)] = m.var('r_cvonm', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:322 — variable r_loan
    r_loan = m.r_loan = IndexedVars()
    for i in s.P:
        for j in s.P2:
            for i0 in s.P:
                for j0 in s.P2:
                    for l in s.L:
                        if _in((i, j),loan_window[i0, j0]):
                            r_loan[(i, j, i0, j0, l)] = m.var('r_loan', (i, j, i0, j0, l), lb=-inf, ub=inf, binary=False)
    rl0ub0d = _maximum(p.r_l0_ub)
    # bm.jl:329 — variable r_l0
    r_l0 = m.r_l0 = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                r_l0[(i, j, l)] = m.var('r_l0', (i, j, l), lb=0.0, ub=inf, binary=False)
    # bm.jl:332 — variable r_l0_d_
    r_l0_d_ = m.r_l0_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        r_l0_d_[(i, j, l, k)] = m.var('r_l0_d_', (i, j, l, k), lb=-inf, ub=inf, binary=False)
    # bm.jl:334 — variable r_l0_pd_
    r_l0_pd_ = m.r_l0_pd_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for _state in (sT, sF):
                    r_l0_pd_[(i, j, l, _state)] = m.var('r_l0_pd_', (i, j, l, _state), lb=0.0, ub=inf, binary=False)
    # bm.jl:337 — variable r_l0add
    r_l0add = m.r_l0add = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                r_l0add[(i, j, l)] = m.var('r_l0add', (i, j, l), lb=0.0, ub=inf, binary=False)
    # bm.jl:343 — variable r_l0add_d_
    r_l0add_d_ = m.r_l0add_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                r_l0add_d_[(i, j, l)] = m.var('r_l0add_d_', (i, j, l), lb=0.0, ub=inf, binary=False)
    # bm.jl:349 — variable r_e_c_d_
    r_e_c_d_ = m.r_e_c_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        r_e_c_d_[(i, j, l, k)] = m.var('r_e_c_d_', (i, j, l, k), lb=0.0, ub=inf, binary=False)
    rleub0d = _maximum(p.r_le_ub)
    # bm.jl:354 — variable r_le
    r_le = m.r_le = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                r_le[(i, j, l)] = m.var('r_le', (i, j, l), lb=0.0, ub=inf, binary=False)
    # bm.jl:359 — variable r_le_d_
    r_le_d_ = m.r_le_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        r_le_d_[(i, j, l, k)] = m.var('r_le_d_', (i, j, l, k), lb=-inf, ub=inf, binary=False)
    # bm.jl:361 — variable r_le_pd_
    r_le_pd_ = m.r_le_pd_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for _state in (sT, sF):
                    r_le_pd_[(i, j, l, _state)] = m.var('r_le_pd_', (i, j, l, _state), lb=0.0, ub=inf, binary=False)
    # bm.jl:366 — variable r_leadd
    r_leadd = m.r_leadd = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                r_leadd[(i, j, l)] = m.var('r_leadd', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:370 — variable r_leadd_d_
    r_leadd_d_ = m.r_leadd_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                r_leadd_d_[(i, j, l)] = m.var('r_leadd_d_', (i, j, l), lb=0.0, ub=inf, binary=False)
    # bm.jl:374 — variable r_le_ped_
    r_le_ped_ = m.r_le_ped_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for _state in (sT, sF):
                    r_le_ped_[(i, j, l, _state)] = m.var('r_le_ped_', (i, j, l, _state), lb=0.0, ub=inf, binary=False)
    # bm.jl:379 — variable r_leadde
    r_leadde = m.r_leadde = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                r_leadde[(i, j, l)] = m.var('r_leadde', (i, j, l), lb=0.0, ub=inf, binary=False)
    # bm.jl:383 — variable r_leadde_d_
    r_leadde_d_ = m.r_leadde_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                r_leadde_d_[(i, j, l)] = m.var('r_leadde_d_', (i, j, l), lb=0.0, ub=inf, binary=False)
    # bm.jl:387 — variable r_pay
    r_pay = m.r_pay = IndexedVars()
    for i in s.P:
        for j in s.P2:
            for l in L:
                r_pay[(i, j, l)] = m.var('r_pay', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:388 — variable r_paye
    r_paye = m.r_paye = IndexedVars()
    for i in s.P:
        for j in s.P2:
            for l in L:
                r_paye[(i, j, l)] = m.var('r_paye', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:402 — variable e_ladd
    e_ladd = m.e_ladd = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                e_ladd[(i, j, l)] = m.var('e_ladd', (i, j, l), lb=0.0, ub=p.e_ladd_ub[l], binary=False)
    # bm.jl:405 — variable e_ladd_d_
    e_ladd_d_ = m.e_ladd_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                e_ladd_d_[(i, j, l)] = m.var('e_ladd_d_', (i, j, l), lb=0.0, ub=p.e_ladd_ub[l], binary=False)
    # bm.jl:408 — variable e_l_d_
    e_l_d_ = m.e_l_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for _state in (sT, sF):
                    e_l_d_[(i, j, l, _state)] = m.var('e_l_d_', (i, j, l, _state), lb=0.0, ub=inf, binary=False)
    # bm.jl:409 — variable e_l_pd_
    e_l_pd_ = m.e_l_pd_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for _state in (sT, sF):
                    e_l_pd_[(i, j, l, _state)] = m.var('e_l_pd_', (i, j, l, _state), lb=0.0, ub=p.e_l_ub[l], binary=False)
    # bm.jl:413 — variable e_ann
    e_ann = m.e_ann = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                e_ann[(i, j, l)] = m.var('e_ann', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:415 — variable e_ann_d_
    e_ann_d_ = m.e_ann_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                e_ann_d_[(i, j, l)] = m.var('e_ann_d_', (i, j, l), lb=0.0, ub=inf, binary=False)
    # bm.jl:416 — variable e_ann_0
    e_ann_0 = m.e_ann_0 = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                e_ann_0[(i, j, l)] = m.var('e_ann_0', (i, j, l), lb=0.0, ub=inf, binary=False)
    # bm.jl:417 — variable e_ann_1
    e_ann_1 = m.e_ann_1 = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                e_ann_1[(i, j, l)] = m.var('e_ann_1', (i, j, l), lb=0.0, ub=inf, binary=False)
    # bm.jl:419 — variable e_pay
    e_pay = m.e_pay = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                e_pay[(i, j, l)] = m.var('e_pay', (i, j, l), lb=0.0, ub=p.e_pay_ub[l], binary=False)
    # bm.jl:422 — variable e_loan
    e_loan = m.e_loan = IndexedVars()
    for i in s.P:
        for j in s.P2:
            for i0 in s.P:
                for j0 in s.P2:
                    for l in s.L:
                        if _in((i, j),loan_window[i0, j0]):
                            e_loan[(i, j, i0, j0, l)] = m.var('e_loan', (i, j, i0, j0, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:426 — variable t_loan_d_
    t_loan_d_ = m.t_loan_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for _state in (sT, sF):
                    t_loan_d_[(i, j, l, _state)] = m.var('t_loan_d_', (i, j, l, _state), lb=0, ub=p.t_loan_bM[l], binary=False)
    # bm.jl:429 — variable t_ret_cost
    t_ret_cost = m.t_ret_cost = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                t_ret_cost[(i, j, l)] = m.var('t_ret_cost', (i, j, l), lb=0.0, ub=p.t_ret_c_bM[l], binary=False)
    # bm.jl:432 — variable t_ret_cost_d_
    t_ret_cost_d_ = m.t_ret_cost_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                t_ret_cost_d_[(i, j, l)] = m.var('t_ret_cost_d_', (i, j, l), lb=0, ub=p.t_ret_c_bM[l], binary=False)
    # bm.jl:438 — variable y_n
    y_n = m.y_n = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        y_n[(i, j, l, k)] = m.var('y_n', (i, j, l, k), lb=-inf, ub=1.0, binary=True)
    # bm.jl:440 — variable n_c0
    n_c0 = m.n_c0 = IndexedVars()
    for i in P:
        for l in L:
            n_c0[(i, l)] = m.var('n_c0', (i, l), lb=0.0, ub=inf, binary=False)
    # bm.jl:442 — variable n_c0_d_
    n_c0_d_ = m.n_c0_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        n_c0_d_[(i, j, l, k)] = m.var('n_c0_d_', (i, j, l, k), lb=0, ub=inf, binary=False)
    # bm.jl:444 — variable n_cp
    n_cp = m.n_cp = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    n_cp[(i, j, l, n)] = m.var('n_cp', (i, j, l, n), lb=0.0, ub=p.n_cp_bM[l, n], binary=False)
    # bm.jl:447 — variable n_cp_d_
    n_cp_d_ = m.n_cp_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if n_filter[l, k]:
                            n_cp_d_[(i, j, l, k, n)] = m.var('n_cp_d_', (i, j, l, k, n), lb=0, ub=inf, binary=False)
    # bm.jl:450 — variable n_l_d_
    n_l_d_ = m.n_l_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        n_l_d_[(i, j, l, k)] = m.var('n_l_d_', (i, j, l, k), lb=-inf, ub=inf, binary=False)
    # bm.jl:451 — variable n_l
    n_l = m.n_l = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                n_l[(i, j, l)] = m.var('n_l', (i, j, l), lb=0, ub=p.n_l_bM[l], binary=False)
    # bm.jl:453 — variable n_ann_d_
    n_ann_d_ = m.n_ann_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        n_ann_d_[(i, j, l, k)] = m.var('n_ann_d_', (i, j, l, k), lb=0, ub=inf, binary=False)
    # bm.jl:455 — variable n_ann
    n_ann = m.n_ann = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                n_ann[(i, j, l)] = m.var('n_ann', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:456 — variable n_ann_0
    n_ann_0 = m.n_ann_0 = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                n_ann_0[(i, j, l)] = m.var('n_ann_0', (i, j, l), lb=0, ub=inf, binary=False)
    # bm.jl:457 — variable n_ann_1
    n_ann_1 = m.n_ann_1 = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                n_ann_1[(i, j, l)] = m.var('n_ann_1', (i, j, l), lb=0, ub=inf, binary=False)
    # bm.jl:460 — variable n_ladd_d_
    n_ladd_d_ = m.n_ladd_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                n_ladd_d_[(i, j, l)] = m.var('n_ladd_d_', (i, j, l), lb=0, ub=p.n_ladd_bM[l], binary=False)
    # bm.jl:461 — variable n_ladd
    n_ladd = m.n_ladd = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                n_ladd[(i, j, l)] = m.var('n_ladd', (i, j, l), lb=0, ub=p.n_l_bM[l], binary=False)
    # bm.jl:463 — variable n_l_pd_
    n_l_pd_ = m.n_l_pd_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for _state in (sT, sF):
                    n_l_pd_[(i, j, l, _state)] = m.var('n_l_pd_', (i, j, l, _state), lb=0, ub=p.n_l_bM[l], binary=False)
    # bm.jl:467 — variable n_loan
    n_loan = m.n_loan = IndexedVars()
    for i in s.P:
        for j in s.P2:
            for i0 in s.P:
                for j0 in s.P2:
                    for l in s.L:
                        if _in((i, j),loan_window[i0, j0]):
                            n_loan[(i, j, i0, j0, l)] = m.var('n_loan', (i, j, i0, j0, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:468 — variable n_pay
    n_pay = m.n_pay = IndexedVars()
    for i in s.P:
        for j in s.P2:
            for l in L:
                n_pay[(i, j, l)] = m.var('n_pay', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:478 — variable n_eh
    n_eh = m.n_eh = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        n_eh[(i, j, l, n)] = m.var('n_eh', (i, j, l, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:479 — variable n_eh_d_
    n_eh_d_ = m.n_eh_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if n_filter[l, k]  and  nd_en_fltr[n]:
                            n_eh_d_[(i, j, l, k, n)] = m.var('n_eh_d_', (i, j, l, k, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:482 — variable n_ehf
    n_ehf = m.n_ehf = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for f in Fu_n[l]:
                    for n in Nd:
                        if nd_en_fltr[n]:
                            n_ehf[(i, j, l, f, n)] = m.var('n_ehf', (i, j, l, f, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:483 — variable n_ehf_d_
    n_ehf_d_ = m.n_ehf_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for f in Fu_n[l]:
                        for n in Nd:
                            if n_filter[l, k]  and  nd_en_fltr[n]:
                                n_ehf_d_[(i, j, l, k, f, n)] = m.var('n_ehf_d_', (i, j, l, k, f, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:486 — variable n_u
    n_u = m.n_u = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        n_u[(i, j, l, n)] = m.var('n_u', (i, j, l, n), lb=0.0, ub=p.n_u_ub[l, n], binary=False)
    # bm.jl:488 — variable n_u_d_
    n_u_d_ = m.n_u_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if n_filter[l, k]  and  nd_en_fltr[n]:
                            n_u_d_[(i, j, l, k, n)] = m.var('n_u_d_', (i, j, l, k, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:492 — variable n_u_onsite
    n_u_onsite = m.n_u_onsite = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        n_u_onsite[(i, j, l, n)] = m.var('n_u_onsite', (i, j, l, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:493 — variable n_u_onsite_d_
    n_u_onsite_d_ = m.n_u_onsite_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if n_filter[l, k]  and  nd_en_fltr[n]:
                            n_u_onsite_d_[(i, j, l, k, n)] = m.var('n_u_onsite_d_', (i, j, l, k, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:497 — variable n_u_ehf
    n_u_ehf = m.n_u_ehf = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for f in Fu_n[l]:
                    for n in Nd:
                        if nd_en_fltr[n]:
                            n_u_ehf[(i, j, l, f, n)] = m.var('n_u_ehf', (i, j, l, f, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:498 — variable n_u_ehf_d_
    n_u_ehf_d_ = m.n_u_ehf_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for f in Fu_n[l]:
                        for n in Nd:
                            if n_filter[l, k]  and  nd_en_fltr[n]:
                                n_u_ehf_d_[(i, j, l, k, f, n)] = m.var('n_u_ehf_d_', (i, j, l, k, f, n), lb=0.0, ub=inf, binary=False)
    # bm.jl:502 — variable n_cpe
    n_cpe = m.n_cpe = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_em_fltr[n]:
                        n_cpe[(i, j, l, n)] = m.var('n_cpe', (i, j, l, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:503 — variable n_cpe_d_
    n_cpe_d_ = m.n_cpe_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if n_filter[l, k]  and  nd_em_fltr[n]:
                            n_cpe_d_[(i, j, l, k, n)] = m.var('n_cpe_d_', (i, j, l, k, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:506 — variable n_fu_e
    n_fu_e = m.n_fu_e = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        n_fu_e[(i, j, l, n)] = m.var('n_fu_e', (i, j, l, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:507 — variable n_fu_e_d_
    n_fu_e_d_ = m.n_fu_e_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if n_filter[l, k]  and  nd_en_fltr[n]:
                            n_fu_e_d_[(i, j, l, k, n)] = m.var('n_fu_e_d_', (i, j, l, k, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:510 — variable n_u_fu_e
    n_u_fu_e = m.n_u_fu_e = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        n_u_fu_e[(i, j, l, n)] = m.var('n_u_fu_e', (i, j, l, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:511 — variable n_u_fu_e_d_
    n_u_fu_e_d_ = m.n_u_fu_e_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if n_filter[l, k]  and  nd_en_fltr[n]:
                            n_u_fu_e_d_[(i, j, l, k, n)] = m.var('n_u_fu_e_d_', (i, j, l, k, n), lb=-inf, ub=inf, binary=False)
    # bm.jl:515 — variable n_ep0
    n_ep0 = m.n_ep0 = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                n_ep0[(i, j, l)] = m.var('n_ep0', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:516 — variable n_ep0_d_
    n_ep0_d_ = m.n_ep0_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        n_ep0_d_[(i, j, l, k)] = m.var('n_ep0_d_', (i, j, l, k), lb=-inf, ub=inf, binary=False)
    # bm.jl:519 — variable n_ep1ge
    n_ep1ge = m.n_ep1ge = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                n_ep1ge[(i, j, l)] = m.var('n_ep1ge', (i, j, l), lb=0.0, ub=p.n_ep1ge_bM[l], binary=False)
    # bm.jl:521 — variable n_ep1ge_d_
    n_ep1ge_d_ = m.n_ep1ge_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        n_ep1ge_d_[(i, j, l, k)] = m.var('n_ep1ge_d_', (i, j, l, k), lb=-inf, ub=inf, binary=False)
    # bm.jl:523 — variable n_ep1gce
    n_ep1gce = m.n_ep1gce = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                n_ep1gce[(i, j, l)] = m.var('n_ep1gce', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:524 — variable n_ep1gce_d_
    n_ep1gce_d_ = m.n_ep1gce_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        n_ep1gce_d_[(i, j, l, k)] = m.var('n_ep1gce_d_', (i, j, l, k), lb=-inf, ub=inf, binary=False)
    # bm.jl:526 — variable n_ep1gcs
    n_ep1gcs = m.n_ep1gcs = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                n_ep1gcs[(i, j, l)] = m.var('n_ep1gcs', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:527 — variable n_ep1gcs_d_
    n_ep1gcs_d_ = m.n_ep1gcs_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        n_ep1gcs_d_[(i, j, l, k)] = m.var('n_ep1gcs_d_', (i, j, l, k), lb=-inf, ub=inf, binary=False)
    # bm.jl:531 — variable n_fstck
    n_fstck = m.n_fstck = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for f in Nf:
                    n_fstck[(i, j, l, f)] = m.var('n_fstck', (i, j, l, f), lb=-inf, ub=inf, binary=False)
    # bm.jl:532 — variable n_fstck_d_
    n_fstck_d_ = m.n_fstck_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for f in Nf:
                        if n_filter[l, k]:
                            n_fstck_d_[(i, j, l, k, f)] = m.var('n_fstck_d_', (i, j, l, k, f), lb=0.0, ub=inf, binary=False)
    # bm.jl:535 — variable n_cfonm_d_
    n_cfonm_d_ = m.n_cfonm_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        n_cfonm_d_[(i, j, l, k)] = m.var('n_cfonm_d_', (i, j, l, k), lb=0.0, ub=inf, binary=False)
    # bm.jl:536 — variable n_cfonm
    n_cfonm = m.n_cfonm = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                n_cfonm[(i, j, l)] = m.var('n_cfonm', (i, j, l), lb=0.0, ub=p.n_cfonm_bM[l], binary=False)
    # bm.jl:539 — variable n_cvonm_d_
    n_cvonm_d_ = m.n_cvonm_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        n_cvonm_d_[(i, j, l, k)] = m.var('n_cvonm_d_', (i, j, l, k), lb=0.0, ub=inf, binary=False)
    # bm.jl:540 — variable n_cvonm
    n_cvonm = m.n_cvonm = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                n_cvonm[(i, j, l)] = m.var('n_cvonm', (i, j, l), lb=0.0, ub=p.n_cvonm_bM[l], binary=False)
    # bm.jl:544 — variable n_u_cost
    n_u_cost = m.n_u_cost = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                n_u_cost[(i, j, l)] = m.var('n_u_cost', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:546 — variable n_ehf_cost
    n_ehf_cost = m.n_ehf_cost = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                n_ehf_cost[(i, j, l)] = m.var('n_ehf_cost', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:552 — variable o_pay
    o_pay = m.o_pay = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                o_pay[(i, j, l)] = m.var('o_pay', (i, j, l), lb=0.0, ub=p.o_pay_bM[l], binary=False)
    # bm.jl:556 — variable o_pay_d_
    o_pay_d_ = m.o_pay_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                o_pay_d_[(i, j, l)] = m.var('o_pay_d_', (i, j, l), lb=0.0, ub=inf, binary=False)
    # bm.jl:558 — variable o_cfonm
    o_cfonm = m.o_cfonm = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                o_cfonm[(i, j, l)] = m.var('o_cfonm', (i, j, l), lb=0.0, ub=p.o_cfonm_bM[l], binary=False)
    # bm.jl:561 — variable o_cfonm_d_
    o_cfonm_d_ = m.o_cfonm_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                o_cfonm_d_[(i, j, l)] = m.var('o_cfonm_d_', (i, j, l), lb=0.0, ub=inf, binary=False)
    # bm.jl:562 — variable o_tcfonm_d_
    o_tcfonm_d_ = m.o_tcfonm_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for _state in (sT, sF):
                    o_tcfonm_d_[(i, j, l, _state)] = m.var('o_tcfonm_d_', (i, j, l, _state), lb=-inf, ub=inf, binary=False)
    # bm.jl:564 — variable o_cvonm
    o_cvonm = m.o_cvonm = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                o_cvonm[(i, j, l)] = m.var('o_cvonm', (i, j, l), lb=0.0, ub=p.o_cvonm_bM[l], binary=False)
    # bm.jl:567 — variable o_cvonm_d_
    o_cvonm_d_ = m.o_cvonm_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                o_cvonm_d_[(i, j, l)] = m.var('o_cvonm_d_', (i, j, l), lb=0.0, ub=inf, binary=False)
    # bm.jl:568 — variable o_tcvonm_d_
    o_tcvonm_d_ = m.o_tcvonm_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for _state in (sT, sF):
                    o_tcvonm_d_[(i, j, l, _state)] = m.var('o_tcvonm_d_', (i, j, l, _state), lb=-inf, ub=inf, binary=False)
    # bm.jl:570 — variable o_cp
    o_cp = m.o_cp = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    o_cp[(i, j, l, n)] = m.var('o_cp', (i, j, l, n), lb=0.0, ub=inf, binary=False)
    # bm.jl:573 — variable o_cp_d_
    o_cp_d_ = m.o_cp_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    o_cp_d_[(i, j, l, n)] = m.var('o_cp_d_', (i, j, l, n), lb=0, ub=inf, binary=False)
    # bm.jl:576 — variable o_tcp_d_
    o_tcp_d_ = m.o_tcp_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for _state in (sT, sF):
                    for n in Nd:
                        o_tcp_d_[(i, j, l, _state, n)] = m.var('o_tcp_d_', (i, j, l, _state, n), lb=0.0, ub=inf, binary=False)
    # bm.jl:577 — variable o_tpay_d_
    o_tpay_d_ = m.o_tpay_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for _state in (sT, sF):
                    o_tpay_d_[(i, j, l, _state)] = m.var('o_tpay_d_', (i, j, l, _state), lb=0.0, ub=inf, binary=False)
    # bm.jl:580 — variable o_rcp_d_
    o_rcp_d_ = m.o_rcp_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for _state in (sT, sF):
                        if r_filter[l, k]:
                            o_rcp_d_[(i, j, l, k, _state)] = m.var('o_rcp_d_', (i, j, l, k, _state), lb=0.0, ub=inf, binary=False)
    # bm.jl:581 — variable o_rcp
    o_rcp = m.o_rcp = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        o_rcp[(i, j, l, k)] = m.var('o_rcp', (i, j, l, k), lb=-inf, ub=inf, binary=False)
    # bm.jl:584 — variable o_u_d_
    o_u_d_ = m.o_u_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        o_u_d_[(i, j, l, n)] = m.var('o_u_d_', (i, j, l, n), lb=0, ub=inf, binary=False)
    # bm.jl:585 — variable o_u
    o_u = m.o_u = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        o_u[(i, j, l, n)] = m.var('o_u', (i, j, l, n), lb=0, ub=p.o_u_ub[l, n], binary=False)
    # bm.jl:587 — variable o_tu_d_
    o_tu_d_ = m.o_tu_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for _state in (sT, sF):
                    for n in Nd:
                        if nd_en_fltr[n]:
                            o_tu_d_[(i, j, l, _state, n)] = m.var('o_tu_d_', (i, j, l, _state, n), lb=0, ub=inf, binary=False)
    # bm.jl:591 — variable o_ehf_d_
    o_ehf_d_ = m.o_ehf_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for f in Fu_r[l]:
                    for n in Nd:
                        if nd_en_fltr[n]:
                            o_ehf_d_[(i, j, l, f, n)] = m.var('o_ehf_d_', (i, j, l, f, n), lb=0, ub=inf, binary=False)
    # bm.jl:592 — variable o_ehf
    o_ehf = m.o_ehf = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for f in Fu_r[l]:
                    for n in Nd:
                        if nd_en_fltr[n]:
                            o_ehf[(i, j, l, f, n)] = m.var('o_ehf', (i, j, l, f, n), lb=0, ub=inf, binary=False)
    # bm.jl:594 — variable o_tehf_d_
    o_tehf_d_ = m.o_tehf_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for f in Fu_r[l]:
                    for _state in (sT, sF):
                        for n in Nd:
                            if nd_en_fltr[n]:
                                o_tehf_d_[(i, j, l, f, _state, n)] = m.var('o_tehf_d_', (i, j, l, f, _state, n), lb=0, ub=inf, binary=False)
    # bm.jl:598 — variable o_ep0_d_
    o_ep0_d_ = m.o_ep0_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                o_ep0_d_[(i, j, l)] = m.var('o_ep0_d_', (i, j, l), lb=0, ub=inf, binary=False)
    # bm.jl:599 — variable o_ep0
    o_ep0 = m.o_ep0 = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                o_ep0[(i, j, l)] = m.var('o_ep0', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:600 — variable o_tep0_d_
    o_tep0_d_ = m.o_tep0_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for _state in (sT, sF):
                    o_tep0_d_[(i, j, l, _state)] = m.var('o_tep0_d_', (i, j, l, _state), lb=0, ub=inf, binary=False)
    # bm.jl:602 — variable o_ep1ge_d_
    o_ep1ge_d_ = m.o_ep1ge_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                o_ep1ge_d_[(i, j, l)] = m.var('o_ep1ge_d_', (i, j, l), lb=0, ub=inf, binary=False)
    # bm.jl:603 — variable o_ep1ge
    o_ep1ge = m.o_ep1ge = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                o_ep1ge[(i, j, l)] = m.var('o_ep1ge', (i, j, l), lb=0.0, ub=p.o_ep1ge_bM[l], binary=False)
    # bm.jl:605 — variable o_tep1ge_d_
    o_tep1ge_d_ = m.o_tep1ge_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for _state in (sT, sF):
                    o_tep1ge_d_[(i, j, l, _state)] = m.var('o_tep1ge_d_', (i, j, l, _state), lb=0, ub=inf, binary=False)
    # bm.jl:608 — variable o_ep1gce_d_
    o_ep1gce_d_ = m.o_ep1gce_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                o_ep1gce_d_[(i, j, l)] = m.var('o_ep1gce_d_', (i, j, l), lb=0, ub=inf, binary=False)
    # bm.jl:609 — variable o_ep1gce
    o_ep1gce = m.o_ep1gce = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                o_ep1gce[(i, j, l)] = m.var('o_ep1gce', (i, j, l), lb=0.0, ub=p.o_ep1gce_bM[l], binary=False)
    # bm.jl:612 — variable o_tep1gce_d_
    o_tep1gce_d_ = m.o_tep1gce_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for _state in (sT, sF):
                    o_tep1gce_d_[(i, j, l, _state)] = m.var('o_tep1gce_d_', (i, j, l, _state), lb=0, ub=inf, binary=False)
    # bm.jl:616 — variable o_ep1gcs_d_
    o_ep1gcs_d_ = m.o_ep1gcs_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                o_ep1gcs_d_[(i, j, l)] = m.var('o_ep1gcs_d_', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:617 — variable o_ep1gcs
    o_ep1gcs = m.o_ep1gcs = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                o_ep1gcs[(i, j, l)] = m.var('o_ep1gcs', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:619 — variable o_tep1gcs_d_
    o_tep1gcs_d_ = m.o_tep1gcs_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for _state in (sT, sF):
                    o_tep1gcs_d_[(i, j, l, _state)] = m.var('o_tep1gcs_d_', (i, j, l, _state), lb=0, ub=inf, binary=False)
    oloanlastub = {l: p.e_loan_ub[l] + _vec(p.r_loan_ub)[l] for l in L}
    # bm.jl:631 — variable o_loan_last
    o_loan_last = m.o_loan_last = IndexedVars()
    for i in P:
        for l in L:
            for _state in (sT, sF):
                o_loan_last[(i, l, _state)] = m.var('o_loan_last', (i, l, _state), lb=0.0, ub=oloanlastub[l], binary=False)
    # bm.jl:635 — variable o_u_cost
    o_u_cost = m.o_u_cost = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                o_u_cost[(i, j, l)] = m.var('o_u_cost', (i, j, l), lb=0.0, ub=inf, binary=False)
    # bm.jl:637 — variable o_ehf_cost
    o_ehf_cost = m.o_ehf_cost = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                o_ehf_cost[(i, j, l)] = m.var('o_ehf_cost', (i, j, l), lb=0.0, ub=inf, binary=False)
    # bm.jl:639 — variable o_ep1gcs_cost
    o_ep1gcs_cost = m.o_ep1gcs_cost = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                o_ep1gcs_cost[(i, j, l)] = m.var('o_ep1gcs_cost', (i, j, l), lb=0.0, ub=inf, binary=False)
    # bm.jl:640 — variable n_ep1gcs_cost
    n_ep1gcs_cost = m.n_ep1gcs_cost = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                n_ep1gcs_cost[(i, j, l)] = m.var('n_ep1gcs_cost', (i, j, l), lb=0.0, ub=inf, binary=False)
    Nd = s.Nd
    Mt = s.Mt
    Ln = s.Ln
    node_mat = p.node_mat
    skip_mb = p.skip_mb
    links_list = p.links_list
    input_mat = p.input_mat
    output_mat = p.output_mat
    # bm.jl:686 — variable r_x_in_d_
    r_x_in_d_ = m.r_x_in_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        for c in Mt:
                            if r_filter[l, k]  and  node_mat[n, c]  and  c  in  input_mat[n]:
                                r_x_in_d_[(i, j, l, k, n, c)] = m.var('r_x_in_d_', (i, j, l, k, n, c), lb=0.0, ub=inf, binary=False)
    # bm.jl:689 — variable r_x_out_d_
    r_x_out_d_ = m.r_x_out_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        for c in Mt:
                            if r_filter[l, k]  and  node_mat[n, c]  and  c  in  output_mat[n]:
                                r_x_out_d_[(i, j, l, k, n, c)] = m.var('r_x_out_d_', (i, j, l, k, n, c), lb=0.0, ub=inf, binary=False)
    # bm.jl:693 — variable r_x_o_s_d_
    r_x_o_s_d_ = m.r_x_o_s_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]:
                            r_x_o_s_d_[(i, j, l, k, n)] = m.var('r_x_o_s_d_', (i, j, l, k, n), lb=0.0, ub=inf, binary=False)
    # bm.jl:696 — variable r_x_in
    r_x_in = m.r_x_in = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  input_mat[n]:
                            r_x_in[(i, j, l, n, c)] = m.var('r_x_in', (i, j, l, n, c), lb=0.0, ub=inf, binary=False)
    # bm.jl:698 — variable r_x_out
    r_x_out = m.r_x_out = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  output_mat[n]:
                            r_x_out[(i, j, l, n, c)] = m.var('r_x_out', (i, j, l, n, c), lb=0.0, ub=inf, binary=False)
    # bm.jl:701 — variable o_x_in
    o_x_in = m.o_x_in = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  input_mat[n]:
                            o_x_in[(i, j, l, n, c)] = m.var('o_x_in', (i, j, l, n, c), lb=0.0, ub=inf, binary=False)
    # bm.jl:703 — variable o_x_in_d_
    o_x_in_d_ = m.o_x_in_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  input_mat[n]:
                            o_x_in_d_[(i, j, l, n, c)] = m.var('o_x_in_d_', (i, j, l, n, c), lb=0.0, ub=inf, binary=False)
    # bm.jl:705 — variable o_tx_in_d_
    o_tx_in_d_ = m.o_tx_in_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        for _state in (sT, sF):
                            if node_mat[n, c]  and  c  in  input_mat[n]:
                                o_tx_in_d_[(i, j, l, n, c, _state)] = m.var('o_tx_in_d_', (i, j, l, n, c, _state), lb=0.0, ub=inf, binary=False)
    # bm.jl:708 — variable o_x_out
    o_x_out = m.o_x_out = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  output_mat[n]:
                            o_x_out[(i, j, l, n, c)] = m.var('o_x_out', (i, j, l, n, c), lb=0.0, ub=inf, binary=False)
    # bm.jl:710 — variable o_x_out_d_
    o_x_out_d_ = m.o_x_out_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  output_mat[n]:
                            o_x_out_d_[(i, j, l, n, c)] = m.var('o_x_out_d_', (i, j, l, n, c), lb=0.0, ub=inf, binary=False)
    # bm.jl:712 — variable o_tx_out_d_
    o_tx_out_d_ = m.o_tx_out_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        for _state in (sT, sF):
                            if node_mat[n, c]  and  c  in  output_mat[n]:
                                o_tx_out_d_[(i, j, l, n, c, _state)] = m.var('o_tx_out_d_', (i, j, l, n, c, _state), lb=0.0, ub=inf, binary=False)
    # bm.jl:717 — variable n_x_in_d_
    n_x_in_d_ = m.n_x_in_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        for c in Mt:
                            if n_filter[l, k]  and  node_mat[n, c]  and  c  in  input_mat[n]:
                                n_x_in_d_[(i, j, l, k, n, c)] = m.var('n_x_in_d_', (i, j, l, k, n, c), lb=0.0, ub=inf, binary=False)
    # bm.jl:720 — variable n_x_out_d_
    n_x_out_d_ = m.n_x_out_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        for c in Mt:
                            if n_filter[l, k]  and  node_mat[n, c]  and  c  in  output_mat[n]:
                                n_x_out_d_[(i, j, l, k, n, c)] = m.var('n_x_out_d_', (i, j, l, k, n, c), lb=0.0, ub=inf, binary=False)
    # bm.jl:725 — variable n_x_o_s_d_
    n_x_o_s_d_ = m.n_x_o_s_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if n_filter[l, k]:
                            n_x_o_s_d_[(i, j, l, k, n)] = m.var('n_x_o_s_d_', (i, j, l, k, n), lb=0.0, ub=inf, binary=False)
    # bm.jl:727 — variable n_x_in
    n_x_in = m.n_x_in = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  input_mat[n]:
                            n_x_in[(i, j, l, n, c)] = m.var('n_x_in', (i, j, l, n, c), lb=0.0, ub=inf, binary=False)
    # bm.jl:729 — variable n_x_out
    n_x_out = m.n_x_out = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  output_mat[n]:
                            n_x_out[(i, j, l, n, c)] = m.var('n_x_out', (i, j, l, n, c), lb=0.0, ub=inf, binary=False)
    # bm.jl:732 — variable r_ups_e_mt_in_d_
    r_ups_e_mt_in_d_ = m.r_ups_e_mt_in_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        r_ups_e_mt_in_d_[(i, j, l, k)] = m.var('r_ups_e_mt_in_d_', (i, j, l, k), lb=0.0, ub=inf, binary=False)
    # bm.jl:733 — variable r_ups_e_mt_in
    r_ups_e_mt_in = m.r_ups_e_mt_in = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                r_ups_e_mt_in[(i, j, l)] = m.var('r_ups_e_mt_in', (i, j, l), lb=0.0, ub=inf, binary=False)
    # bm.jl:735 — variable o_ups_e_mt_in_d_
    o_ups_e_mt_in_d_ = m.o_ups_e_mt_in_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                o_ups_e_mt_in_d_[(i, j, l)] = m.var('o_ups_e_mt_in_d_', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:736 — variable o_ups_e_mt_in
    o_ups_e_mt_in = m.o_ups_e_mt_in = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                o_ups_e_mt_in[(i, j, l)] = m.var('o_ups_e_mt_in', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:737 — variable o_tups_e_mt_in_d_
    o_tups_e_mt_in_d_ = m.o_tups_e_mt_in_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for _state in (sT, sF):
                    o_tups_e_mt_in_d_[(i, j, l, _state)] = m.var('o_tups_e_mt_in_d_', (i, j, l, _state), lb=0, ub=inf, binary=False)
    # bm.jl:739 — variable n_ups_e_mt_in_d_
    n_ups_e_mt_in_d_ = m.n_ups_e_mt_in_d_ = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        n_ups_e_mt_in_d_[(i, j, l, k)] = m.var('n_ups_e_mt_in_d_', (i, j, l, k), lb=0.0, ub=inf, binary=False)
    # bm.jl:740 — variable n_ups_e_mt_in
    n_ups_e_mt_in = m.n_ups_e_mt_in = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                n_ups_e_mt_in[(i, j, l)] = m.var('n_ups_e_mt_in', (i, j, l), lb=0.0, ub=inf, binary=False)
    # bm.jl:742 — variable o_x_in_cost
    o_x_in_cost = m.o_x_in_cost = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                o_x_in_cost[(i, j, l)] = m.var('o_x_in_cost', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:743 — variable n_x_in_cost
    n_x_in_cost = m.n_x_in_cost = IndexedVars()
    for i in P:
        for j in P2:
            for l in L:
                n_x_in_cost[(i, j, l)] = m.var('n_x_in_cost', (i, j, l), lb=-inf, ub=inf, binary=False)
    # bm.jl:750 — constraint o_logic_1_y_i_
    for i in P:
        for j in P2:
            for l in L:
                if j<n_subperiods:
                    m.add((y_o[i, j+1, l]) - (y_o[i, j, l]), '<=', 0, name='o_logic_1_y_i_' + ":" + repr((i, j, l)))
    fP = min(s.P)
    fP2 = min(P2)
    fKr = min(Kr)
    fKn = min(Kn)
    # bm.jl:770 — constraint r_logic_budget_s
    for i in P:
        for j in P2:
            for l in L:
                if (i,j) != (fP,fP2):
                    m.add((sum(y_r[i, j, l, k] for k in Kr if (k > fKr  and  r_filter[l, k]))) - (1), '<=', 0, name='r_logic_budget_s' + ":" + repr((i, j, l)))
    # bm.jl:776 — constraint r_logic_tier_1_0m_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if (k>fKr  and  r_filter[l, k]):
                        m.add((1) - (y_r[i, j, l, k] + y_r[i, j, l, fKr]), '>=', 0, name='r_logic_tier_1_0m_e_' + ":" + repr((i, j, l, k)))
    # bm.jl:781 — constraint r_logic_budget_y_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if k>fKr  and  j<n_subperiods  and  r_filter[l, k]:
                        m.add((y_r[i, j+1, l, k]) - (y_r[i, j, l, k]), '>=', 0, name='r_logic_budget_y_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:786 — constraint r_logic_onoff_1_y_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if j<n_subperiods  and  r_filter[l, k]:
                        m.add((y_o[i, j, l] + 1 - y_r[i, j, l, k] + y_r[i, j+1, l, k]) - (1), '>=', 0, name='r_logic_onoff_1_y_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:791 — constraint r_logic_onoff_2_y_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if j<n_subperiods  and  r_filter[l, k]:
                        m.add((y_o[i, j, l] + 1 - y_r[i, j+1, l, k] + y_r[i, j, l, k]) - (1), '>=', 0, name='r_logic_onoff_2_y_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:806 — constraint exp_d0_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((e_c_d_[i, j, l, sT]) - (p.e_C[l] * x_d_[i, j, l, sT]), '==', 0, name='exp_d0_e_' + ":" + repr((i, j, l)))
    # bm.jl:810 — constraint exp_ncw_m_i0_
    for i in P:
        for j in P2:
            for l in L:
                m.add((e_c_d_[i, j, l, sT]) - (p.e_c_ub[l] * y_e[i, j, l]), '<=', 0, name='exp_ncw_m_i0_' + ":" + repr((i, j, l)))
    # bm.jl:813 — constraint exp_ncw_s_
    for i in P:
        for j in P2:
            for l in L:
                m.add((e_c[i, j, l]) - (e_c_d_[i, j, l, sT]), '==', 0, name='exp_ncw_s_' + ":" + repr((i, j, l)))
    # bm.jl:817 — constraint exp_x_m_i0_
    for i in P:
        for j in P2:
            for l in L:
                m.add((x_d_[i, j, l, sT]) - (p.x_ub * y_e[i, j, l]), '<=', 0, name='exp_x_m_i0_' + ":" + repr((i, j, l)))
    # bm.jl:820 — constraint exp_x_m_i3_
    for i in P:
        for j in P2:
            for l in L:
                m.add((x_d_[i, j, l, sT]) - (1.0 * y_e[i, j, l]), '>=', 0, name='exp_x_m_i3_' + ":" + repr((i, j, l)))
    # bm.jl:824 — constraint exp_x_m_i1_
    for i in P:
        for j in P2:
            for l in L:
                m.add((x_d_[i, j, l, sF]) - (p.x_ub * (1 - y_e[i, j, l])), '<=', 0, name='exp_x_m_i1_' + ":" + repr((i, j, l)))
    # bm.jl:827 — constraint exp_x_s_
    for i in P:
        for j in P2:
            for l in L:
                m.add((x[i, l]) - (x_d_[i, j, l, sT] + x_d_[i, j, l, sF]), '==', 0, name='exp_x_s_' + ":" + repr((i, j, l)))
    # bm.jl:833 — constraint e_l_d0_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((e_l_d_[i, j, l, sT]) - (p.e_loanFact[l] * e_c_d_[i, j, l, sT]), '==', 0, name='e_l_d0_e_' + ":" + repr((i, j, l)))
    # bm.jl:842 — constraint e_l_m_i0_
    for i in P:
        for j in P2:
            for l in L:
                m.add((e_l_d_[i, j, l, sT]) - (p.e_l_ub[l] * y_e[i, j, l]), '<=', 0, name='e_l_m_i0_' + ":" + repr((i, j, l)))
    # bm.jl:850 — constraint e_l_s_
    for i in P:
        for j in P2:
            for l in L:
                m.add((e_l[i, j, l]) - (e_l_d_[i, j, l, sT]), '==', 0, name='e_l_s_' + ":" + repr((i, j, l)))
    # bm.jl:856 — constraint e_ann_d0_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((e_ann_d_[i, j, l]) - (p.e_Ann[l] * e_l_d_[i, j, l, sT]), '==', 0, name='e_ann_d0_e_' + ":" + repr((i, j, l)))
    # bm.jl:860 — constraint e_ann_m_i0_
    for i in P:
        for j in P2:
            for l in L:
                m.add((e_ann_d_[i, j, l]) - (p.e_ann_ub[l] * y_e[i, j, l]), '<=', 0, name='e_ann_m_i0_' + ":" + repr((i, j, l)))
    # bm.jl:863 — constraint e_ann_s_
    for i in P:
        for j in P2:
            for l in L:
                m.add((e_ann[i, j, l]) - (e_ann_d_[i, j, l]), '==', 0, name='e_ann_s_' + ":" + repr((i, j, l)))
    # bm.jl:871 — constraint e_logic_1_y_i_
    for i in P:
        for j in P2:
            for l in L:
                if j<n_subperiods:
                    m.add((y_e[i, j+1, l]) - (y_e[i, j, l]), '>=', 0, name='e_logic_1_y_i_' + ":" + repr((i, j, l)))
    # bm.jl:876 — constraint e_logic_onoff_exp_1_i_
    for i in P:
        for j in P2:
            for l in L:
                if j<n_subperiods:
                    m.add((y_o[i, j, l] + 1 - y_e[i, j, l] + y_e[i, j+1, l]) - (1), '>=', 0, name='e_logic_onoff_exp_1_i_' + ":" + repr((i, j, l)))
    # bm.jl:880 — constraint e_logic_onoff_exp_2_i_
    for i in P:
        for j in P2:
            for l in L:
                if j<n_subperiods:
                    m.add((y_o[i, j, l] + 1 - y_e[i, j+1, l] + y_e[i, j, l]) - (1), '>=', 0, name='e_logic_onoff_exp_2_i_' + ":" + repr((i, j, l)))
    # bm.jl:891 — constraint e_ladd_d0_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((e_ladd_d_[i, j, l]) - (e_l_pd_[i, j, l, sT]), '==', 0, name='e_ladd_d0_e_' + ":" + repr((i, j, l)))
    # bm.jl:894 — constraint e_ladd_m_0_y_i_
    for i in P:
        for j in P2:
            for l in L:
                if j>fP2:
                    m.add((e_ladd_d_[i, j, l]) - (p.e_ladd_ub[l] * (y_e[i, j, l] - y_e[i, j-1, l])), '<=', 0, name='e_ladd_m_0_y_i_' + ":" + repr((i, j, l)))
    # bm.jl:900 — constraint e_l_m_1_y_i_
    for i in P:
        for j in P2:
            for l in L:
                if j>fP2:
                    m.add((e_l_pd_[i, j, l, sF]) - (p.e_l_ub[l] * (1 - y_e[i, j, l] + y_e[i, j-1, l])), '<=', 0, name='e_l_m_1_y_i_' + ":" + repr((i, j, l)))
    # bm.jl:907 — constraint e_ladd_s_e
    for i in P:
        for j in P2:
            for l in L:
                m.add((e_ladd[i, j, l]) - (e_ladd_d_[i, j, l]), '==', 0, name='e_ladd_s_e' + ":" + repr((i, j, l)))
    # bm.jl:911 — constraint e_l_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((e_l[i, j, l]) - (e_l_pd_[i, j, l, sT] + e_l_pd_[i, j, l, sF]), '==', 0, name='e_l_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:920 — constraint def_e_pay
    for i0 in s.P:
        for j0 in s.P2:
            for l in L:
                m.add((e_pay[i0, j0, l]) - (e_ladd[i0, j0, l] / loan_subperiods), '==', 0, name='def_e_pay' + ":" + repr((i0, j0, l)))
    # bm.jl:926 — constraint e_loan_bal_intrap_0
    for i in s.P:
        for j in s.P2:
            for i0 in s.P:
                for j0 in s.P2:
                    for l in L:
                        if _in((i, j),loan_window[i0, j0])  and  _in((i, j), pay_window[i0, j0])  and  j<n_subperiods  and  (i0,j0)==(i,j):
                            m.add((e_loan[i, j+1, i0, j0, l]) - (- e_pay[i0, j0, l] + e_ladd[i0, j0, l]), '==', 0, name='e_loan_bal_intrap_0' + ":" + repr((i, j, i0, j0, l)))
    # bm.jl:936 — constraint e_loan_bal_intrap
    for i in s.P:
        for j in s.P2:
            for i0 in s.P:
                for j0 in s.P2:
                    for l in L:
                        if _in((i, j),loan_window[i0, j0])  and  _in((i, j), pay_window[i0, j0])  and  j<n_subperiods  and  (i0,j0)!=(i,j):
                            m.add((e_loan[i, j+1, i0, j0, l]) - (e_loan[i, j, i0, j0, l] - e_pay[i0, j0, l]), '==', 0, name='e_loan_bal_intrap' + ":" + repr((i, j, i0, j0, l)))
    # bm.jl:946 — constraint e_loan_bal_intra_unpay
    for i in s.P:
        for j in s.P2:
            for i0 in s.P:
                for j0 in s.P2:
                    for l in L:
                        if _in((i,j), loan_window[i0, j0])  and  not ((i, j) in pay_window[i0, j0])  and  j<n_subperiods:
                            m.add((e_loan[i, j+1, i0, j0, l]) - (e_loan[i, j, i0, j0, l]), '==', 0, name='e_loan_bal_intra_unpay' + ":" + repr((i, j, i0, j0, l)))
    # bm.jl:955 — constraint e_loan_bal_per_0
    for i in s.P:
        for i0 in s.P:
            for j0 in s.P2:
                for l in L:
                    if _in((i, n_subperiods),loan_window[i0, j0])  and  _in((i, n_subperiods), pay_window[i0, j0])  and  i < n_periods  and  (i0, j0) == (i,n_subperiods):
                        m.add((e_loan[i+1, 1, i0, j0, l]) - (- e_pay[i0, j0, l] + e_ladd[i0, j0, l]), '==', 0, name='e_loan_bal_per_0' + ":" + repr((i, i0, j0, l)))
    # bm.jl:966 — constraint e_loan_bal_per
    for i in s.P:
        for i0 in s.P:
            for j0 in s.P2:
                for l in L:
                    if _in((i, n_subperiods),loan_window[i0, j0])  and  _in((i, n_subperiods), pay_window[i0, j0])  and  i < n_periods  and  (i0,j0) != (i,n_subperiods):
                        m.add((e_loan[i+1, 1, i0, j0, l]) - (e_loan[i, n_subperiods, i0, j0, l] - e_pay[i0, j0, l]), '==', 0, name='e_loan_bal_per' + ":" + repr((i, i0, j0, l)))
    # bm.jl:977 — constraint e_loan_bal_per_unpay
    for i in s.P:
        for i0 in s.P:
            for j0 in s.P2:
                for l in L:
                    if _in((i,n_subperiods), loan_window[i0, j0])  and  not _in((i,n_subperiods), pay_window[i0, j0])  and  i<n_periods:
                        m.add((e_loan[i+1, 1, i0, j0, l]) - (e_loan[i, n_subperiods, i0, j0, l]), '==', 0, name='e_loan_bal_per_unpay' + ":" + repr((i, i0, j0, l)))
    # bm.jl:990 — constraint cpb_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((cpb[i, j, l]) - (p.c0[l] + e_c[i, j, l]), '==', 0, name='cpb_e_' + ":" + repr((i, j, l)))
    # bm.jl:1000 — constraint r_cp_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]:
                            m.add((r_cp_d_[i, j, l, k, n]) - (r_cpb_d_[i, j, l, k] * p.r_Kkey_j[n, k]), '<=', 0, name='r_cp_d_e_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:1005 — constraint r_cp_ub_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        m.add((r_cpb_d_[i, j, l, k] * p.min_cpr[l]) - (r_cp_d_[i, j, l, k, key_node]), '<=', 0, name='r_cp_ub_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:1010 — constraint r_cp_lb_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]:
                            m.add((r_cp_d_[i, j, l, k, n]) - (p.r_cp_ub[l, n] * y_r[i, j, l, k]), '<=', 0, name='r_cp_lb_i_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:1013 — constraint r_cp_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    m.add((r_cp[i, j, l, n]) - (sum(r_cp_d_[i, j, l, k, n] for k in Kr if r_filter[l, k])), '==', 0, name='r_cp_s_e_' + ":" + repr((i, j, l, n)))
    # bm.jl:1017 — constraint r_cp_x_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]:
                            m.add((r_cp_d_[i, j, l, k, n]) - (r_x_out_d_[i, j, l, k, n, ckey[n]] + r_x_o_s_d_[i, j, l, k, n]), '==', 0, name='r_cp_x_e_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:1023 — constraint r_cpb_d_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        m.add((r_cpb_d_[i, j, l, k]) - (p.r_cpb_ub[l] * y_r[i, j, l, k]), '<=', 0, name='r_cpb_d_m_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:1026 — constraint r_cpb_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((cpb[i, j, l]) - (sum(r_cpb_d_[i, j, l, k] for k in Kr if r_filter[l, k])), '==', 0, name='r_cpb_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1032 — constraint r_eh_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]  and  nd_en_fltr[n]:
                            m.add((r_eh_d_[i, j, l, k, n]) - (p.r_c_H[l, k, n] * r_cp_d_[i, j, l, k, n] + p.r_rhs_H[l, k, n] * y_r[i, j, l, k]), '==', 0, name='r_eh_d_e_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:1038 — constraint r_eh_d_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]  and  nd_en_fltr[n]:
                            m.add((r_eh_d_[i, j, l, k, n]) - (p.r_eh_ub[l, n] * y_r[i, j, l, k]), '<=', 0, name='r_eh_d_m_i_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:1042 — constraint r_eh_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        m.add((r_eh[i, j, l, n]) - (sum(r_eh_d_[i, j, l, k, n] for k in Kr if r_filter[l, k])), '==', 0, name='r_eh_s_e_' + ":" + repr((i, j, l, n)))
    # bm.jl:1049 — constraint r_ehf_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for f in Fu_r[l]:
                        for n in Nd:
                            if r_filter[l, k]  and  nd_en_fltr[n]:
                                m.add((r_ehf_d_[i, j, l, k, f, n]) - (p.r_c_F[f, l, k, n] * r_eh_d_[i, j, l, k, n] + p.r_rhs_F[f, l, k, n] * y_r[i, j, l, k]), '==', 0, name='r_ehf_d_e_' + ":" + repr((i, j, l, k, f, n)))
    # bm.jl:1056 — constraint r_ehf_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for f in Fu_r[l]:
                        for n in Nd:
                            if r_filter[l, k]  and  nd_en_fltr[n]:
                                m.add((r_ehf_d_[i, j, l, k, f, n]) - (p.r_ehf_ub[f, l, n] * y_r[i, j, l, k]), '<=', 0, name='r_ehf_m_i_' + ":" + repr((i, j, l, k, f, n)))
    # bm.jl:1060 — constraint r_ehf_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for f in Fu_r[l]:
                    for n in Nd:
                        if nd_en_fltr[n]:
                            m.add((r_ehf[i, j, l, f, n]) - (sum(r_ehf_d_[i, j, l, k, f, n] for k in Kr if r_filter[l, k])), '==', 0, name='r_ehf_s_e_' + ":" + repr((i, j, l, f, n)))
    # bm.jl:1067 — constraint r_u_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]  and  nd_en_fltr[n]:
                            m.add((r_u_d_[i, j, l, k, n]) - ((1-p.r_c_UonSite[i, j, l, k, n]) * (p.r_c_U[l, k, n] * r_cp_d_[i, j, l, k, n] + p.r_rhs_U[l, k, n] * y_r[i,j,l,k])), '==', 0, name='r_u_d_e_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:1074 — constraint r_u_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]  and  nd_en_fltr[n]:
                            m.add((r_u_d_[i, j, l, k, n]) - (p.r_u_ub[l, n] * y_r[i, j, l, k]), '<=', 0, name='r_u_i_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:1078 — constraint r_u_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        m.add((r_u[i, j, l, n]) - (sum(r_u_d_[i, j, l, k, n] for k in Kr if r_filter[l, k])), '==', 0, name='r_u_s_e_' + ":" + repr((i, j, l, n)))
    # bm.jl:1084 — constraint r_u_onsite_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]  and  nd_en_fltr[n]:
                            m.add((r_u_onsite_d_[i, j, l, k, n]) - ((p.r_c_UonSite[i, j, l, k, n]) * (p.r_c_U[l, k, n] * r_cp_d_[i, j, l, k, n] + p.r_rhs_U[l,k,n] * y_r[i,j,l,k])), '==', 0, name='r_u_onsite_d_e_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:1090 — constraint r_u_onsite_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]  and  nd_en_fltr[n]:
                            m.add((r_u_onsite_d_[i, j, l, k, n]) - (p.r_u_ub[l, n] * y_r[i, j, l, k]), '<=', 0, name='r_u_onsite_i_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:1094 — constraint r_u_onsite_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        m.add((r_u_onsite[i, j, l, n]) - (sum(r_u_onsite_d_[i, j, l, k, n] for k in Kr if r_filter[l, k])), '==', 0, name='r_u_onsite_s_e_' + ":" + repr((i, j, l, n)))
    # bm.jl:1100 — constraint r_u_ehf_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for f in Fu_r[l]:
                        for n in Nd:
                            if r_filter[l, k]  and  nd_en_fltr[n]:
                                m.add((r_u_ehf_d_[i, j, l, k, f, n]) - (p.r_c_Hr[f, l, k, n] * p.r_c_Fgenf[f, l, k, n] * r_u_onsite_d_[i, j, l, k, n]), '==', 0, name='r_u_ehf_e_' + ":" + repr((i, j, l, k, f, n)))
    # bm.jl:1106 — constraint r_u_ehf_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for f in Fu_r[l]:
                        for n in Nd:
                            if r_filter[l, k]  and  nd_en_fltr[n]:
                                m.add((r_u_ehf_d_[i, j, l, k, f, n]) - (p.r_u_ehf_ub[f, l, n] * y_r[i, j, l, k]), '<=', 0, name='r_u_ehf_i_' + ":" + repr((i, j, l, k, f, n)))
    # bm.jl:1111 — constraint r_u_ehf_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for f in Fu_r[l]:
                        for n in Nd:
                            if r_filter[l, k]  and  nd_en_fltr[n]:
                                m.add((r_u_ehf[i, j, l, f, n]) - (sum(r_u_ehf_d_[i, j, l, k, f, n] for k in Kr if r_filter[l, k])), '==', 0, name='r_u_ehf_s_e_' + ":" + repr((i, j, l, k, f, n)))
    # bm.jl:1120 — constraint r_cpe_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]  and  nd_em_fltr[n]:
                            m.add((r_cpe_d_[i, j, l, k, n]) - (p.r_c_cpe[l, k, n] * r_cp_d_[i, j, l, k, n] + p.r_rhs_cpe[l, k, n] * y_r[i, j, l, k]), '==', 0, name='r_cpe_d_e_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:1127 — constraint r_cpe_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]  and  nd_em_fltr[n]:
                            m.add((r_cpe_d_[i, j, l, k, n]) - (p.r_cpe_ub[l, n] * y_r[i, j, l, k]), '<=', 0, name='r_cpe_m_i_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:1131 — constraint r_cpe_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_em_fltr[n]:
                        m.add((r_cpe[i, j, l, n]) - (sum(r_cpe_d_[i, j, l, k, n] for k in Kr if r_filter[l, k])), '==', 0, name='r_cpe_s_e_' + ":" + repr((i, j, l, n)))
    # bm.jl:1137 — constraint r_fu_e_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]  and  nd_en_fltr[n]:
                            m.add((r_fu_e_d_[i, j, l, k, n]) - (sum(p.r_c_Fe[f, l, k, n] * r_ehf_d_[i, j, l, k, f, n] for f in Fu_r[l])), '==', 0, name='r_fu_e_d_e_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:1143 — constraint r_fu_e_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]  and  nd_en_fltr[n]:
                            m.add((r_fu_e_d_[i, j, l, k, n]) - (p.r_fu_e_ub[l, n] * y_r[i, j, l, k]), '<=', 0, name='r_fu_e_m_i_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:1148 — constraint r_fu_e_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        m.add((r_fu_e[i, j, l, n]) - (sum(r_fu_e_d_[i, j, l, k, n] for k in Kr if r_filter[l, k])), '==', 0, name='r_fu_e_s_e_' + ":" + repr((i, j, l, n)))
    # bm.jl:1153 — constraint r_u_fu_e_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]  and  nd_en_fltr[n]:
                            m.add((r_u_fu_e_d_[i, j, l, k, n]) - (sum(p.r_c_Fe[f, l, k, n] * r_u_ehf_d_[i, j, l, k, f, n] for f in Fu_r[l])), '==', 0, name='r_u_fu_e_d_e_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:1159 — constraint r_u_fu_e_d_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        if r_filter[l, k]  and  nd_en_fltr[n]:
                            m.add((r_u_fu_e_d_[i, j, l, k, n]) - (p.r_u_fu_e_ub[l, n] * y_r[i, j, l, k]), '<=', 0, name='r_u_fu_e_d_m_i_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:1164 — constraint r_u_fu_e_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        m.add((r_u_fu_e[i, j, l, n]) - (sum(r_u_fu_e_d_[i, j, l, k, n] for k in Kr if r_filter[l, k])), '==', 0, name='r_u_fu_e_s_e_' + ":" + repr((i, j, l, n)))
    # bm.jl:1179 — constraint r_ep0_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        m.add((r_ep0_d_[i, j, l, k]) - (sum(r_fu_e_d_[i, j, l, k, n] for n in Nd if nd_en_fltr[n]) + sum(r_cpe_d_[i, j, l, k, n] for n in Nd if nd_em_fltr[n]) + sum(r_u_fu_e_d_[i, j, l, k, n] for n in Nd if nd_en_fltr[n])), '==', 0, name='r_ep0_d_e_' + ":" + repr((i, j, l, k)))
    # bm.jl:1189 — constraint r_ep0_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        m.add((r_ep0_d_[i, j, l, k]) - (p.r_ep0_ub[l] * y_r[i, j, l, k]), '<=', 0, name='r_ep0_m_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:1192 — constraint r_ep0_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_ep0[i, j, l]) - (sum(r_ep0_d_[i, j, l, k] for k in Kr if r_filter[l, k])), '==', 0, name='r_ep0_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1199 — constraint r_ep1ge_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        m.add((r_ep1ge_d_[i, j, l, k]) - ((1.0 - p.r_chi[l, k]) * r_ep0_d_[i, j, l, k]), '==', 0, name='r_ep1ge_d_e_' + ":" + repr((i, j, l, k)))
    # bm.jl:1204 — constraint r_ep1ge_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        m.add((r_ep1ge_d_[i, j, l, k]) - (p.r_ep1ge_ub[l] * y_r[i, j, l, k]), '<=', 0, name='r_ep1ge_m_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:1207 — constraint r_ep1ge_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_ep1ge[i, j, l]) - (sum(r_ep1ge_d_[i, j, l, k] for k in Kr if r_filter[l, k])), '==', 0, name='r_ep1ge_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1214 — constraint r_ep1gce_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        m.add((r_ep1gce_d_[i, j, l, k]) - (p.r_chi[l, k] * (1 - p.r_sigma[l, k]) * r_ep0_d_[i, j, l, k]), '==', 0, name='r_ep1gce_d_e_' + ":" + repr((i, j, l, k)))
    # bm.jl:1218 — constraint r_ep1gce_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        m.add((r_ep1gce_d_[i, j, l, k]) - (p.r_ep1gce_ub[l] * y_r[i, j, l, k]), '<=', 0, name='r_ep1gce_m_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:1221 — constraint r_ep1gce_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_ep1gce[i, j, l]) - (sum(r_ep1gce_d_[i, j, l, k] for k in Kr if r_filter[l, k])), '==', 0, name='r_ep1gce_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1228 — constraint r_ep1gcs_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        m.add((r_ep1gcs_d_[i, j, l, k]) - (p.r_chi[l, k] * p.r_sigma[l, k] * r_ep0_d_[i, j, l, k]), '==', 0, name='r_ep1gcs_d_e_' + ":" + repr((i, j, l, k)))
    # bm.jl:1233 — constraint r_ep1gcs_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        m.add((r_ep1gcs_d_[i, j, l, k]) - (p.r_ep1gcs_ub[l] * y_r[i, j, l, k]), '<=', 0, name='r_ep1gcs_m_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:1237 — constraint r_ep1gcs_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_ep1gcs[i, j, l]) - (sum(r_ep1gcs_d_[i, j, l, k] for k in Kr if r_filter[l, k])), '==', 0, name='r_ep1gcs_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1243 — constraint r_ups_e_mt_in_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        m.add((r_ups_e_mt_in_d_[i, j, l, k]) - (sum( sum(p.r_c_upsein_rate[l, c] * r_x_in_d_[i, j, l, k, n, c] for c in Mt if node_mat[n, c]  and  c  in  input_mat[n] ) for n in Nd)), '==', 0, name='r_ups_e_mt_in_e_' + ":" + repr((i, j, l, k)))
    # bm.jl:1252 — constraint r_ups_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        m.add((r_ups_e_mt_in_d_[i, j, l, k]) - (p.r_ups_e_mt_in_ub[l] * y_r[i, j, l, k]), '<=', 0, name='r_ups_m_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:1255 — constraint r_ups_e_mt_in_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_ups_e_mt_in[i, j, l]) - (sum(r_ups_e_mt_in_d_[i, j, l, k] for k in Kr if r_filter[l, k])), '==', 0, name='r_ups_e_mt_in_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1260 — constraint o_tups_e_mt_in_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_ups_e_mt_in[i, j, l]) - (o_tups_e_mt_in_d_[i, j, l, sT] + o_tups_e_mt_in_d_[i, j, l, sF]), '==', 0, name='o_tups_e_mt_in_e_' + ":" + repr((i, j, l)))
    # bm.jl:1266 — constraint o_tups_e_mt_in_m1_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_tups_e_mt_in_d_[i, j, l, sT]) - (p.o_ups_e_mt_in_ub[l] * y_o[i, j, l]), '<=', 0, name='o_tups_e_mt_in_m1_i_' + ":" + repr((i, j, l)))
    # bm.jl:1271 — constraint o_tups_e_mt_in_m0_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_tups_e_mt_in_d_[i, j, l, sF]) - (p.o_ups_e_mt_in_ub[l] * (1 - y_o[i, j, l])), '<=', 0, name='o_tups_e_mt_in_m0_i_' + ":" + repr((i, j, l)))
    # bm.jl:1276 — constraint o_ups_e_mt_d_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_ups_e_mt_in_d_[i, j, l]) - (o_tups_e_mt_in_d_[i, j, l, sT]), '==', 0, name='o_ups_e_mt_d_e_' + ":" + repr((i, j, l)))
    # bm.jl:1279 — constraint o_ups_e_mt_d_m_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_ups_e_mt_in_d_[i, j, l]) - (p.o_ups_e_mt_in_ub[l] * y_o[i, j, l]), '<=', 0, name='o_ups_e_mt_d_m_i_' + ":" + repr((i, j, l)))
    # bm.jl:1282 — constraint o_ups_e_mt_d_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_ups_e_mt_in[i, j, l]) - (o_ups_e_mt_in_d_[i, j, l]), '==', 0, name='o_ups_e_mt_d_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1320 — constraint r_cfonm_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        m.add((r_cfonm_d_[i, j, l, k]) - (p.r_c_fOnm[l, k] * r_cpb_d_[i, j, l, k] + p.r_rhs_fOnm[l, k] * y_r[i, j, l, k]), '==', 0, name='r_cfonm_d_e_' + ":" + repr((i, j, l, k)))
    # bm.jl:1325 — constraint r_cfonm_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        m.add((r_cfonm_d_[i, j, l, k]) - (p.r_cfonm_ub[l] * y_r[i, j, l, k]), '<=', 0, name='r_cfonm_m_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:1328 — constraint r_cfonm_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_cfonm[i, j, l]) - (sum(r_cfonm_d_[i, j, l, k] for k in Kr if r_filter[l, k])), '==', 0, name='r_cfonm_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1334 — constraint r_cvonm_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        m.add((r_cvonm_d_[i, j, l, k]) - (p.r_c_vOnm[l, k] * r_cp_d_[i, j, l, k, key_node] + p.r_rhs_vOnm[l, k] * y_r[i, j, l, k]), '==', 0, name='r_cvonm_d_e_' + ":" + repr((i, j, l, k)))
    # bm.jl:1339 — constraint r_cvonm_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        m.add((r_cvonm_d_[i, j, l, k]) - (p.r_cvonm_ub[l] * y_r[i, j, l, k]), '<=', 0, name='r_cvonm_m_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:1342 — constraint r_cvonm_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_cvonm[i, j, l]) - (sum(r_cvonm_d_[i, j, l, k] for k in Kr if r_filter[l, k])), '==', 0, name='r_cvonm_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1351 — constraint r_ec_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        m.add((e_c[i, j, l]) - (sum(r_e_c_d_[i, j, l, k] for k in Kr if r_filter[l, k])), '==', 0, name='r_ec_s_e_' + ":" + repr((i, j, l, k)))
    # bm.jl:1356 — constraint r_ec_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        m.add((r_e_c_d_[i, j, l, k]) - (p.r_e_c_ub[l] * y_r[i, j, l, k]), '<=', 0, name='r_ec_m_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:1362 — constraint r_l0_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        m.add((r_l0_d_[i, j, l, k]) - (p.r_loanFact[i, j, l, k] * p.c0[l] * y_r[i, j, l, k]), '==', 0, name='r_l0_d_e_' + ":" + repr((i, j, l, k)))
    # bm.jl:1371 — constraint r_l0_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_l0[i, j, l]) - (sum(r_l0_d_[i, j, l, k] for k in Kr if r_filter[l, k])), '==', 0, name='r_l0_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1377 — constraint r_le_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        m.add((r_le_d_[i, j, l, k]) - (p.r_loanFact[i, j, l, k] * r_e_c_d_[i, j, l, k]), '==', 0, name='r_le_d_e_' + ":" + repr((i, j, l, k)))
    # bm.jl:1382 — constraint r_le_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    if r_filter[l, k]:
                        m.add((r_le_d_[i, j, l, k]) - (p.r_le_ub[l, k] * y_r[i, j, l, k]), '<=', 0, name='r_le_m_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:1386 — constraint r_le_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_le[i, j, l]) - (sum(r_le_d_[i, j, l, k] for k in Kr if r_filter[l, k])), '==', 0, name='r_le_s_e_' + ":" + repr((i, j, l)))
    rl0ub1d = _maximum(p.r_l0_ub, dims=2)
    # bm.jl:1427 — constraint r_l0add_d_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_l0add_d_[i, j, l]) - (r_l0_pd_[i, j, l, sT]), '==', 0, name='r_l0add_d_e_' + ":" + repr((i, j, l)))
    # bm.jl:1431 — constraint r_l0add_m_0_y_i_
    for i in P:
        for j in P2:
            for l in L:
                if j>fP2:
                    m.add((r_l0add_d_[i, j, l]) - (rl0ub1d[l] * (y_r[i, j-1, l, fKr] - y_r[i, j, l, fKr])), '<=', 0, name='r_l0add_m_0_y_i_' + ":" + repr((i, j, l)))
    # bm.jl:1437 — constraint r_l0add_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_l0add[i, j, l]) - (r_l0add_d_[i, j, l]), '==', 0, name='r_l0add_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1442 — constraint r_l0_m_1_y_i_
    for i in P:
        for j in P2:
            for l in L:
                if j>fP2:
                    m.add((r_l0_pd_[i, j, l, sF]) - (rl0ub1d[l] * (1.0 - y_r[i, j-1, l, fKr] + y_r[i, j, l, fKr])), '<=', 0, name='r_l0_m_1_y_i_' + ":" + repr((i, j, l)))
    # bm.jl:1447 — constraint r_l0_p_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_l0[i, j, l]) - (r_l0_pd_[i, j, l, sT] + r_l0_pd_[i, j, l, sF]), '==', 0, name='r_l0_p_s_e_' + ":" + repr((i, j, l)))
    rleub1d = _maximum(p.r_le_ub, dims=2)
    # bm.jl:1455 — constraint r_leadd_d_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_leadd_d_[i, j, l]) - (r_le_pd_[i, j, l, sT]), '==', 0, name='r_leadd_d_e_' + ":" + repr((i, j, l)))
    # bm.jl:1459 — constraint r_leadd_m_1_y_i_
    for i in P:
        for j in P2:
            for l in L:
                if j>fP2:
                    m.add((r_leadd_d_[i, j, l]) - (rleub1d[l] * (y_r[i, j-1, l, fKr] - y_r[i, j, l, fKr])), '<=', 0, name='r_leadd_m_1_y_i_' + ":" + repr((i, j, l)))
    # bm.jl:1464 — constraint r_leadd_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_leadd[i, j, l]) - (r_leadd_d_[i, j, l]), '==', 0, name='r_leadd_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1468 — constraint r_le_m_1_y_i_
    for i in P:
        for j in P2:
            for l in L:
                if j>fP2:
                    m.add((r_le_pd_[i, j, l, sF]) - (rleub1d[l] * (1.0 - y_r[i, j-1, l, fKr] + y_r[i, j, l, fKr])), '<=', 0, name='r_le_m_1_y_i_' + ":" + repr((i, j, l)))
    # bm.jl:1473 — constraint r_le_p_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_le[i, j, l]) - (r_le_pd_[i, j, l, sT] + r_le_pd_[i, j, l, sF]), '==', 0, name='r_le_p_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1478 — constraint r_leadde_d_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_leadde_d_[i, j, l]) - (r_le_ped_[i, j, l, sT]), '==', 0, name='r_leadde_d_e_' + ":" + repr((i, j, l)))
    # bm.jl:1481 — constraint r_leadde_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_leadde[i, j, l]) - (r_leadde_d_[i, j, l]), '==', 0, name='r_leadde_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1484 — constraint r_le_pe_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_le[i, j, l]) - (r_le_ped_[i, j, l, sT] + r_le_ped_[i, j, l, sF]), '==', 0, name='r_le_pe_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1487 — constraint r_leadd_m_2_y_i_
    for i in P:
        for j in P2:
            for l in L:
                if j>fP2:
                    m.add((r_leadde_d_[i, j, l]) - (rleub1d[l] * (y_e[i, j, l] - y_e[i, j-1, l])), '<=', 0, name='r_leadd_m_2_y_i_' + ":" + repr((i, j, l)))
    # bm.jl:1491 — constraint r_le_m_2_y_i_
    for i in P:
        for j in P2:
            for l in L:
                if j>fP2:
                    m.add((r_le_ped_[i, j, l, sF]) - (rleub1d[l] * (1 - y_e[i, j, l] + y_e[i, j-1, l])), '<=', 0, name='r_le_m_2_y_i_' + ":" + repr((i, j, l)))
    # bm.jl:1508 — constraint def_r_pay
    for i0 in s.P:
        for j0 in s.P2:
            for l in L:
                m.add((r_pay[i0, j0, l]) - (r_l0add[i0, j0, l] / loan_subperiods), '==', 0, name='def_r_pay' + ":" + repr((i0, j0, l)))
    # bm.jl:1511 — constraint def_r_paye_
    for i0 in s.P:
        for j0 in s.P2:
            for l in L:
                m.add((r_paye[i0, j0, l]) - (r_leadd[i0, j0, l] / loan_subperiods), '==', 0, name='def_r_paye_' + ":" + repr((i0, j0, l)))
    # bm.jl:1517 — constraint loan_bal_intrap_0
    for i in s.P:
        for j in s.P2:
            for i0 in s.P:
                for j0 in s.P2:
                    for l in L:
                        if _in((i, j),loan_window[i0, j0])  and  _in((i, j), pay_window[i0, j0])  and  j<n_subperiods  and  (i0,j0)==(i,j):
                            m.add((r_loan[i, j+1, i0, j0, l]) - (- r_pay[i0, j0, l] - r_paye[i0, j0, l] + r_l0add[i0, j0, l] + r_leadd[i0, j0, l]), '==', 0, name='loan_bal_intrap_0' + ":" + repr((i, j, i0, j0, l)))
    # bm.jl:1529 — constraint loan_bal_intrap
    for i in s.P:
        for j in s.P2:
            for i0 in s.P:
                for j0 in s.P2:
                    for l in L:
                        if _in((i, j),loan_window[i0, j0])  and  _in((i, j), pay_window[i0, j0])  and  j<n_subperiods  and  (i0,j0)!=(i,j):
                            m.add((r_loan[i, j+1, i0, j0, l]) - (r_loan[i, j, i0, j0, l] - r_pay[i0, j0, l] - r_paye[i0, j0, l]), '==', 0, name='loan_bal_intrap' + ":" + repr((i, j, i0, j0, l)))
    # bm.jl:1540 — constraint loan_bal_intra_unpay
    for i in s.P:
        for j in s.P2:
            for i0 in s.P:
                for j0 in s.P2:
                    for l in L:
                        if _in((i,j), loan_window[i0, j0])  and  not ((i, j) in pay_window[i0, j0])  and  j<n_subperiods:
                            m.add((r_loan[i, j+1, i0, j0, l]) - (r_loan[i, j, i0, j0, l]), '==', 0, name='loan_bal_intra_unpay' + ":" + repr((i, j, i0, j0, l)))
    # bm.jl:1549 — constraint loan_bal_per_0
    for i in s.P:
        for i0 in s.P:
            for j0 in s.P2:
                for l in L:
                    if _in((i, n_subperiods),loan_window[i0, j0])  and  _in((i, n_subperiods), pay_window[i0, j0])  and  i < n_periods  and  (i0, j0) == (i,n_subperiods):
                        m.add((r_loan[i+1, 1, i0, j0, l]) - (- r_pay[i0, j0, l] - r_paye[i0, j0, l] + r_l0add[i0, j0, l] + r_leadd[i0, j0, l]), '==', 0, name='loan_bal_per_0' + ":" + repr((i, i0, j0, l)))
    # bm.jl:1562 — constraint loan_bal_per
    for i in s.P:
        for i0 in s.P:
            for j0 in s.P2:
                for l in L:
                    if _in((i, n_subperiods),loan_window[i0, j0])  and  _in((i, n_subperiods), pay_window[i0, j0])  and  i < n_periods  and  (i0,j0) != (i,n_subperiods):
                        m.add((r_loan[i+1, 1, i0, j0, l]) - (r_loan[i, n_subperiods, i0, j0, l] - r_pay[i0, j0, l] - r_paye[i0, j0, l]), '==', 0, name='loan_bal_per' + ":" + repr((i, i0, j0, l)))
    # bm.jl:1574 — constraint loan_bal_per_unpay
    for i in s.P:
        for i0 in s.P:
            for j0 in s.P2:
                for l in L:
                    if _in((i,n_subperiods), loan_window[i0, j0])  and  not _in((i,n_subperiods), pay_window[i0, j0])  and  i<n_periods:
                        m.add((r_loan[i+1, 1, i0, j0, l]) - (r_loan[i, n_subperiods, i0, j0, l]), '==', 0, name='loan_bal_per_unpay' + ":" + repr((i, i0, j0, l)))
    # bm.jl:1587 — constraint r_mass_b_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        for c in Mt:
                            if r_filter[l, k]  and  node_mat[n, c]  and  not skip_mb[n, c]:
                                m.add((f_xin(m, i, j, l, k, n, c, input_mat, mode="r") - f_xout(m, i, j, l, k, n, c, output_mat, mode="r") + (p.r_Kmb[n, c, ckey[n], k] * r_cp_d_[i, j, l, k, n])) - (0.0), '==', 0, name='r_mass_b_e_' + ":" + repr((i, j, l, k, n, c)))
    # bm.jl:1597 — constraint r_comp_link_e_
    for i in P:
        for j in P2:
            for l in L:
                for lnk in Ln:
                    m.add((r_x_out[i, j, l, links_list[lnk][1], links_list[lnk][3]]) - (r_x_in[i, j, l, links_list[lnk][2], links_list[lnk][3]]), '==', 0, name='r_comp_link_e_' + ":" + repr((i, j, l, lnk)))
    # bm.jl:1603 — constraint r_x_in_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        for c in Mt:
                            if r_filter[l, k]  and  node_mat[n, c]  and  c  in  input_mat[n]:
                                m.add((r_x_in_d_[i, j, l, k, n, c]) - (p.r_x_in_ub[l, n, c] * y_r[i, j, l, k]), '<=', 0, name='r_x_in_m_i_' + ":" + repr((i, j, l, k, n, c)))
    # bm.jl:1609 — constraint r_x_in_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  input_mat[n]:
                            m.add((r_x_in[i, j, l, n, c]) - (sum(r_x_in_d_[i, j, l, k, n, c] for k in Kr if r_filter[l, k])), '==', 0, name='r_x_in_s_e_' + ":" + repr((i, j, l, n, c)))
    # bm.jl:1615 — constraint r_x_out_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kr:
                    for n in Nd:
                        for c in Mt:
                            if r_filter[l, k]  and  node_mat[n, c]  and  c  in  output_mat[n]:
                                m.add((r_x_out_d_[i, j, l, k, n, c]) - (p.r_x_out_ub[l, n, c] * y_r[i, j, l, k]), '<=', 0, name='r_x_out_m_i_' + ":" + repr((i, j, l, k, n, c)))
    # bm.jl:1621 — constraint r_x_out_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  output_mat[n]:
                            m.add((r_x_out[i, j, l, n, c]) - (sum(r_x_out_d_[i, j, l, k, n, c] for k in Kr if r_filter[l, k])), '==', 0, name='r_x_out_s_e_' + ":" + repr((i, j, l, n, c)))
    # bm.jl:1628 — constraint o_tx_in_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  input_mat[n]:
                            m.add((r_x_in[i, j, l, n, c]) - (o_tx_in_d_[i, j, l, n, c, sT] + o_tx_in_d_[i, j, l, n, c, sF]), '==', 0, name='o_tx_in_s_e_' + ":" + repr((i, j, l, n, c)))
    # bm.jl:1633 — constraint o_tx_in_m1_i_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  input_mat[n]:
                            m.add((o_tx_in_d_[i, j, l, n, c, sT]) - (p.o_x_in_ub[l, n, c] * y_o[i, j, l]), '<=', 0, name='o_tx_in_m1_i_' + ":" + repr((i, j, l, n, c)))
    # bm.jl:1638 — constraint o_tx_in_m0_i_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  input_mat[n]:
                            m.add((o_tx_in_d_[i, j, l, n, c, sF]) - (p.o_x_in_ub[l, n, c] * (1 - y_o[i, j, l])), '<=', 0, name='o_tx_in_m0_i_' + ":" + repr((i, j, l, n, c)))
    # bm.jl:1643 — constraint o_x_in_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  input_mat[n]:
                            m.add((o_x_in_d_[i, j, l, n, c]) - (o_tx_in_d_[i, j, l, n, c, sT]), '==', 0, name='o_x_in_d_e_' + ":" + repr((i, j, l, n, c)))
    # bm.jl:1647 — constraint o_x_in_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  input_mat[n]:
                            m.add((o_x_in_d_[i, j, l, n, c]) - (p.o_x_in_ub[l, n, c] * y_o[i, j, l]), '<=', 0, name='o_x_in_m_i_' + ":" + repr((i, j, l, n, c)))
    # bm.jl:1652 — constraint o_x_in_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  input_mat[n]:
                            m.add((o_x_in[i, j, l, n, c]) - (o_x_in_d_[i, j, l, n, c]), '==', 0, name='o_x_in_s_e_' + ":" + repr((i, j, l, n, c)))
    # bm.jl:1658 — constraint o_tx_out_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  output_mat[n]:
                            m.add((r_x_out[i, j, l, n, c]) - (o_tx_out_d_[i, j, l, n, c, sT] + o_tx_out_d_[i, j, l, n, c, sF]), '==', 0, name='o_tx_out_s_e_' + ":" + repr((i, j, l, n, c)))
    # bm.jl:1664 — constraint o_tx_out_m1_i_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  output_mat[n]:
                            m.add((o_tx_out_d_[i, j, l, n, c, sT]) - (p.o_x_out_ub[l, n, c] * y_o[i, j, l]), '<=', 0, name='o_tx_out_m1_i_' + ":" + repr((i, j, l, n, c)))
    # bm.jl:1670 — constraint o_tx_out_m0_i_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  output_mat[n]:
                            m.add((o_tx_out_d_[i, j, l, n, c, sF]) - (p.o_x_out_ub[l, n, c] * (1 - y_o[i, j, l])), '<=', 0, name='o_tx_out_m0_i_' + ":" + repr((i, j, l, n, c)))
    # bm.jl:1676 — constraint o_x_out_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  output_mat[n]:
                            m.add((o_x_out_d_[i, j, l, n, c]) - (o_tx_out_d_[i, j, l, n, c, sT]), '==', 0, name='o_x_out_d_e_' + ":" + repr((i, j, l, n, c)))
    # bm.jl:1681 — constraint o_x_out_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  output_mat[n]:
                            m.add((o_x_out_d_[i, j, l, n, c]) - (p.o_x_out_ub[l, n, c] * y_o[i, j, l]), '<=', 0, name='o_x_out_m_i_' + ":" + repr((i, j, l, n, c)))
    # bm.jl:1687 — constraint o_x_out_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  output_mat[n]:
                            m.add((o_x_out[i, j, l, n, c]) - (o_x_out_d_[i, j, l, n, c]), '==', 0, name='o_x_out_s_e_' + ":" + repr((i, j, l, n, c)))
    # bm.jl:1700 — constraint t_ret_c_d_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((t_ret_cost_d_[i, j, l]) - (t_loan_d_[i, j, l, sT]), '==', 0, name='t_ret_c_d_e_' + ":" + repr((i, j, l)))
    # bm.jl:1705 — constraint t_ret_c_bm_0_y_i_
    for i in P:
        for j in P2:
            for l in L:
                if j<n_subperiods:
                    m.add((t_ret_cost_d_[i, j, l]) - (p.t_ret_c_bM[l] * (y_o[i, j, l] - y_o[i, j+1, l])), '<=', 0, name='t_ret_c_bm_0_y_i_' + ":" + repr((i, j, l)))
    # bm.jl:1713 — constraint r_loan_d_bm_0_y_i_
    for i in P:
        for j in P2:
            for l in L:
                if j<n_subperiods:
                    m.add((t_loan_d_[i, j, l, sT]) - (p.t_loan_bM[l] * (y_o[i, j, l] - y_o[i,j+1, l])), '<=', 0, name='r_loan_d_bm_0_y_i_' + ":" + repr((i, j, l)))
    # bm.jl:1719 — constraint r_loan_d_bm_1_y_i_
    for i in P:
        for j in P2:
            for l in L:
                if j<n_subperiods:
                    m.add((t_loan_d_[i, j, l, sF]) - (p.t_loan_bM[l] * (1 + y_o[i,j+1, l] - y_o[i, j, l])), '<=', 0, name='r_loan_d_bm_1_y_i_' + ":" + repr((i, j, l)))
    # bm.jl:1724 — constraint t_ret_c_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((t_ret_cost[i, j, l]) - (t_ret_cost_d_[i, j, l]), '==', 0, name='t_ret_c_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1729 — constraint r_loan_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((sum(r_loan[i, j, i0, j0, l] for j0 in P2 for i0 in P if _in((i, j),loan_window[i0, j0])  and  _ne((i,j),(i0,j0)) ) + sum(e_loan[i, j, i0, j0, l] for j0 in P2 for i0 in P if _in((i, j),loan_window[i0, j0])  and  _ne((i,j),(i0,j0)) ) + r_l0add[i, j, l] + r_leadd[i, j, l] + e_ladd[i, j, l]) - (t_loan_d_[i, j, l, sT] + t_loan_d_[i, j, l, sF]), '==', 0, name='r_loan_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1742 — constraint o_pay_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_pay[i, j, l]) - (o_pay_d_[i, j, l]), '==', 0, name='o_pay_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1745 — constraint o_pay_d1_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_pay_d_[i, j, l]) - (o_tpay_d_[i, j, l, sT]), '==', 0, name='o_pay_d1_e_' + ":" + repr((i, j, l)))
    # bm.jl:1749 — constraint o_pay_m_1_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_pay_d_[i, j, l]) - (p.o_pay_bM[l] * y_o[i, j, l]), '<=', 0, name='o_pay_m_1_i_' + ":" + repr((i, j, l)))
    # bm.jl:1753 — constraint o_tpay_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((sum(r_pay[i0, j0, l] for i0 in P for j0 in P2 if _in((i,j),pay_window[i0, j0])) + sum(r_paye[i0, j0, l] for i0 in P for j0 in P2 if _in((i,j),pay_window[i0, j0])) +sum(e_pay[i0, j0, l] for i0 in P for j0 in P2 if _in((i,j),pay_window[i0, j0]))) - (o_tpay_d_[i, j, l, sF] + o_tpay_d_[i, j, l, sT]), '==', 0, name='o_tpay_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1760 — constraint o_tpay_m1_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_tpay_d_[i, j, l, sT]) - (p.o_pay_bM[l] * y_o[i, j, l]), '<=', 0, name='o_tpay_m1_i_' + ":" + repr((i, j, l)))
    # bm.jl:1764 — constraint o_tpay_m0_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_tpay_d_[i, j, l, sF]) - (p.o_pay_bM[l] * (1 - y_o[i, j, l])), '<=', 0, name='o_tpay_m0_i_' + ":" + repr((i, j, l)))
    # bm.jl:1768 — constraint o_cfonm_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_cfonm[i, j, l]) - (o_cfonm_d_[i, j, l]), '==', 0, name='o_cfonm_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1771 — constraint o_cfonm_d1_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_cfonm_d_[i, j, l]) - (o_tcfonm_d_[i, j, l, sT]), '==', 0, name='o_cfonm_d1_e_' + ":" + repr((i, j, l)))
    # bm.jl:1775 — constraint o_cfonm_m_1_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_cfonm_d_[i, j, l]) - (p.o_cfonm_bM[l] * y_o[i, j, l]), '<=', 0, name='o_cfonm_m_1_i_' + ":" + repr((i, j, l)))
    # bm.jl:1779 — constraint o_tcfonm_m1_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_tcfonm_d_[i, j, l, sT]) - (p.o_cfonm_bM[l] * y_o[i, j, l]), '<=', 0, name='o_tcfonm_m1_i_' + ":" + repr((i, j, l)))
    # bm.jl:1782 — constraint o_tcfonm_m0_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_tcfonm_d_[i, j, l, sF]) - (p.o_cfonm_bM[l] * (1 - y_o[i, j, l])), '<=', 0, name='o_tcfonm_m0_i_' + ":" + repr((i, j, l)))
    # bm.jl:1785 — constraint o_tcfonm_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_cfonm[i, j, l]) - (o_tcfonm_d_[i, j, l, sT] + o_tcfonm_d_[i, j, l, sF]), '==', 0, name='o_tcfonm_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1790 — constraint o_cvonm_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_cvonm[i, j, l]) - (o_cvonm_d_[i, j, l]), '==', 0, name='o_cvonm_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1793 — constraint o_cvonm_d1_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_cvonm_d_[i, j, l]) - (o_tcvonm_d_[i, j, l, sT]), '==', 0, name='o_cvonm_d1_e_' + ":" + repr((i, j, l)))
    # bm.jl:1797 — constraint o_cvonm_m_1_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_cvonm_d_[i, j, l]) - (p.o_cvonm_bM[l] * y_o[i, j, l]), '<=', 0, name='o_cvonm_m_1_i_' + ":" + repr((i, j, l)))
    # bm.jl:1801 — constraint o_tcvonm_m1_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_tcvonm_d_[i, j, l, sT]) - (p.o_cvonm_bM[l] * y_o[i, j, l]), '<=', 0, name='o_tcvonm_m1_i_' + ":" + repr((i, j, l)))
    # bm.jl:1804 — constraint o_tcvonm_m0_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_tcvonm_d_[i, j, l, sF]) - (p.o_cvonm_bM[l] * (1 - y_o[i, j, l])), '<=', 0, name='o_tcvonm_m0_i_' + ":" + repr((i, j, l)))
    # bm.jl:1807 — constraint o_tcvonm_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_cvonm[i, j, l]) - (o_tcvonm_d_[i, j, l, sT] + o_tcvonm_d_[i, j, l, sF]), '==', 0, name='o_tcvonm_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1819 — constraint n_cp0_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    m.add((n_cp_d_[i, j, l, fKn, n]) - (0), '==', 0, name='n_cp0_d_e_' + ":" + repr((i, j, l, n)))
    # bm.jl:1822 — constraint n_cp_ub_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if k>fKn  and  n_filter[l, k]:
                            m.add((n_cp_d_[i, j, l, k, n]) - (n_c0_d_[i, j, l, k] * p.n_Kkey_j[n, k]), '<=', 0, name='n_cp_ub_i_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:1827 — constraint n_cp_lb_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if k>fKn  and  n_filter[l, k]:
                            m.add((n_c0_d_[i, j, l, k] * p.min_cpr[l]) - (n_cp_d_[i, j, l, k, key_node]), '<=', 0, name='n_cp_lb_i_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:1837 — constraint n_cp_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if k>fKn  and  n_filter[l, k]:
                            m.add((n_cp_d_[i, j, l, k, n]) - (p.n_cp_bM[l, n] * y_n[i, j, l, k]), '<=', 0, name='n_cp_m_i_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:1842 — constraint n_cap_add_d_lo_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        m.add((n_c0_d_[i, j, l, k]) - (p.n_c0_lo[l, k] * y_n[i, j, l, k]), '>=', 0, name='n_cap_add_d_lo_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:1846 — constraint n_cap_add_d_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        m.add((n_c0_d_[i, j, l, k]) - (p.n_c0_bM[l] * y_n[i, j, l, k]), '<=', 0, name='n_cap_add_d_m_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:1851 — constraint n_cp_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    m.add((n_cp[i, j, l, n]) - (sum(n_cp_d_[i, j, l, k, n] for k in Kn if n_filter[l, k])), '==', 0, name='n_cp_s_e_' + ":" + repr((i, j, l, n)))
    # bm.jl:1856 — constraint n_cp_x_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if n_filter[l, k]:
                            m.add((n_cp_d_[i, j, l, k, n]) - (n_x_out_d_[i, j, l, k, n, ckey[n]] + n_x_o_s_d_[i, j, l, k, n]), '==', 0, name='n_cp_x_e_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:1861 — constraint n_c0_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((n_c0[i, l]) - (sum(n_c0_d_[i, j, l, k] for k in Kn if n_filter[l, k])), '==', 0, name='n_c0_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1871 — constraint n_l_d0_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        m.add((n_l_d_[i, j, l, k]) - (p.n_loanFact[l, k] * n_c0_d_[i, j, l, k]), '==', 0, name='n_l_d0_e_' + ":" + repr((i, j, l, k)))
    # bm.jl:1876 — constraint n_l_m_i0_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        m.add((n_l_d_[i, j, l, k]) - (p.n_l_bM[l] * y_n[i, j, l, k]), '<=', 0, name='n_l_m_i0_' + ":" + repr((i, j, l, k)))
    # bm.jl:1879 — constraint n_l_s_
    for i in P:
        for j in P2:
            for l in L:
                m.add((n_l[i, j, l]) - (sum(n_l_d_[i, j, l, k] for k in Kn if n_filter[l, k])), '==', 0, name='n_l_s_' + ":" + repr((i, j, l)))
    # bm.jl:1892 — constraint n_ladd_m_0_y_i_
    for i in P:
        for j in P2:
            for l in L:
                if j>fP2:
                    m.add((n_ladd_d_[i, j, l]) - (p.n_ladd_bM[l] * (y_n[i, j-1, l, fKn] - y_n[i, j, l, fKn])), '<=', 0, name='n_ladd_m_0_y_i_' + ":" + repr((i, j, l)))
    # bm.jl:1897 — constraint n_l_pd_m_1_y_i_
    for i in P:
        for j in P2:
            for l in L:
                if j>fP2:
                    m.add((n_l_pd_[i, j, l, sF]) - (p.n_l_bM[l] * (1 + y_n[i, j, l, fKn] - y_n[i, j-1, l, fKn])), '<=', 0, name='n_l_pd_m_1_y_i_' + ":" + repr((i, j, l)))
    # bm.jl:1902 — constraint n_ladd_d0_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((n_ladd_d_[i, j, l]) - (n_l_pd_[i, j, l, sT]), '==', 0, name='n_ladd_d0_e_' + ":" + repr((i, j, l)))
    # bm.jl:1905 — constraint n_ladd_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((n_ladd[i, j, l]) - (n_ladd_d_[i, j, l]), '==', 0, name='n_ladd_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1908 — constraint n_l_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((n_l[i, j, l]) - (n_l_pd_[i, j, l, sT] + n_l_pd_[i, j, l, sF]), '==', 0, name='n_l_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:1914 — constraint def_n_pay
    for i0 in s.P:
        for j0 in s.P2:
            for l in L:
                m.add((n_pay[i0, j0, l]) - (n_ladd[i0, j0, l] / loan_subperiods), '==', 0, name='def_n_pay' + ":" + repr((i0, j0, l)))
    # bm.jl:1920 — constraint n_loan_bal_intrap_0
    for i in s.P:
        for j in s.P2:
            for i0 in s.P:
                for j0 in s.P2:
                    for l in L:
                        if _in((i, j),loan_window[i0, j0])  and  _in((i, j), pay_window[i0, j0])  and  j<n_subperiods  and  (i0,j0)==(i,j):
                            m.add((n_loan[i, j+1, i0, j0, l]) - (- n_pay[i0, j0, l] + n_ladd[i0, j0, l]), '==', 0, name='n_loan_bal_intrap_0' + ":" + repr((i, j, i0, j0, l)))
    # bm.jl:1930 — constraint n_loan_bal_intrap
    for i in s.P:
        for j in s.P2:
            for i0 in s.P:
                for j0 in s.P2:
                    for l in L:
                        if _in((i, j),loan_window[i0, j0])  and  _in((i, j), pay_window[i0, j0])  and  j<n_subperiods  and  (i0,j0)!=(i,j):
                            m.add((n_loan[i, j+1, i0, j0, l]) - (n_loan[i, j, i0, j0, l] - n_pay[i0, j0, l]), '==', 0, name='n_loan_bal_intrap' + ":" + repr((i, j, i0, j0, l)))
    # bm.jl:1940 — constraint n_loan_bal_intra_unpay
    for i in s.P:
        for j in s.P2:
            for i0 in s.P:
                for j0 in s.P2:
                    for l in L:
                        if _in((i,j), loan_window[i0, j0])  and  not ((i, j) in pay_window[i0, j0])  and  j<n_subperiods:
                            m.add((n_loan[i, j+1, i0, j0, l]) - (n_loan[i, j, i0, j0, l]), '==', 0, name='n_loan_bal_intra_unpay' + ":" + repr((i, j, i0, j0, l)))
    # bm.jl:1949 — constraint n_loan_bal_per_0
    for i in s.P:
        for i0 in s.P:
            for j0 in s.P2:
                for l in L:
                    if _in((i, n_subperiods),loan_window[i0, j0])  and  _in((i, n_subperiods), pay_window[i0, j0])  and  i < n_periods  and  (i0, j0) == (i,n_subperiods):
                        m.add((n_loan[i+1, 1, i0, j0, l]) - (- n_pay[i0, j0, l] + n_ladd[i0, j0, l]), '==', 0, name='n_loan_bal_per_0' + ":" + repr((i, i0, j0, l)))
    # bm.jl:1960 — constraint n_loan_bal_per
    for i in s.P:
        for i0 in s.P:
            for j0 in s.P2:
                for l in L:
                    if _in((i, n_subperiods),loan_window[i0, j0])  and  _in((i, n_subperiods), pay_window[i0, j0])  and  i < n_periods  and  (i0,j0) != (i,n_subperiods):
                        m.add((n_loan[i+1, 1, i0, j0, l]) - (n_loan[i, n_subperiods, i0, j0, l] - n_pay[i0, j0, l]), '==', 0, name='n_loan_bal_per' + ":" + repr((i, i0, j0, l)))
    # bm.jl:1971 — constraint n_loan_bal_per_unpay
    for i in s.P:
        for i0 in s.P:
            for j0 in s.P2:
                for l in L:
                    if _in((i,n_subperiods), loan_window[i0, j0])  and  not _in((i,n_subperiods), pay_window[i0, j0])  and  i<n_periods:
                        m.add((n_loan[i+1, 1, i0, j0, l]) - (n_loan[i, n_subperiods, i0, j0, l]), '==', 0, name='n_loan_bal_per_unpay' + ":" + repr((i, i0, j0, l)))
    # bm.jl:1986 — constraint n_eh_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if n_filter[l, k]  and  nd_en_fltr[n]:
                            m.add((n_eh_d_[i, j, l, k, n]) - (p.n_c_H[l, k, n] * n_cp_d_[i, j, l, k, n] + p.n_rhs_H[l, k, n] * y_n[i, j, l, k]), '==', 0, name='n_eh_d_e_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:1993 — constraint n_eh_d_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if n_filter[l, k]  and  nd_en_fltr[n]:
                            m.add((n_eh_d_[i, j, l, k, n]) - (p.n_eh_ub[l, n] * y_n[i, j, l, k]), '<=', 0, name='n_eh_d_m_i_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:1998 — constraint n_eh_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        m.add((n_eh[i, j, l, n]) - (sum(n_eh_d_[i, j, l, k, n] for k in Kn if n_filter[l, k])), '==', 0, name='n_eh_s_e_' + ":" + repr((i, j, l, n)))
    # bm.jl:2004 — constraint n_ehf_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for f in Fu_n[l]:
                        for n in Nd:
                            if n_filter[l, k]  and  nd_en_fltr[n]:
                                m.add((n_ehf_d_[i, j, l, k, f, n]) - (p.n_c_F[f, l, k, n] * n_eh_d_[i, j, l, k, n] + p.n_rhs_F[f, l, k, n] * y_n[i, j, l, k]), '==', 0, name='n_ehf_d_e_' + ":" + repr((i, j, l, k, f, n)))
    # bm.jl:2012 — constraint n_ehf_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for f in Fu_n[l]:
                        for n in Nd:
                            if n_filter[l, k]  and  nd_en_fltr[n]:
                                m.add((n_ehf_d_[i, j, l, k, f, n]) - (p.n_ehf_ub[f, l, n] * y_n[i, j, l, k]), '<=', 0, name='n_ehf_m_i_' + ":" + repr((i, j, l, k, f, n)))
    # bm.jl:2017 — constraint n_ehf_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for f in Fu_n[l]:
                    for n in Nd:
                        if nd_en_fltr[n]:
                            m.add((n_ehf[i, j, l, f, n]) - (sum(n_ehf_d_[i, j, l, k, f, n] for k in Kn if n_filter[l, k])), '==', 0, name='n_ehf_s_e_' + ":" + repr((i, j, l, f, n)))
    # bm.jl:2025 — constraint n_u_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if n_filter[l, k]  and  nd_en_fltr[n]:
                            m.add((n_u_d_[i, j, l, k, n]) - ((1-p.n_c_UonSite[i, j, l, k, n]) * (p.n_c_U[l, k, n] * n_cp_d_[i, j, l, k, n] + p.n_rhs_U[l, k, n] * y_n[i, j, l, k])), '==', 0, name='n_u_d_e_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:2032 — constraint n_u_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if n_filter[l, k]  and  nd_en_fltr[n]:
                            m.add((n_u_d_[i, j, l, k, n]) - (p.n_u_ub[l, n] * y_n[i, j, l, k]), '<=', 0, name='n_u_i_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:2037 — constraint n_u_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        m.add((n_u[i, j, l, n]) - (sum(n_u_d_[i, j, l, k, n] for k in Kn if n_filter[l, k])), '==', 0, name='n_u_s_e_' + ":" + repr((i, j, l, n)))
    # bm.jl:2043 — constraint n_u_onsite_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if n_filter[l, k]  and  nd_en_fltr[n]:
                            m.add((n_u_onsite_d_[i, j, l, k, n]) - ((p.n_c_UonSite[i, j, l, k, n]) * (p.n_c_U[l, k, n] * n_cp_d_[i, j, l, k, n] + p.n_rhs_U[l, k, n] * y_n[i, j, l, k])), '==', 0, name='n_u_onsite_d_e_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:2050 — constraint n_u_onsite_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if n_filter[l, k]  and  nd_en_fltr[n]:
                            m.add((n_u_onsite_d_[i, j, l, k, n]) - (p.n_u_ub[l, n] * y_n[i, j, l, k]), '<=', 0, name='n_u_onsite_i_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:2055 — constraint n_u_onsite_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        m.add((n_u_onsite[i, j, l, n]) - (sum(n_u_onsite_d_[i, j, l, k, n] for k in Kn if n_filter[l, k])), '==', 0, name='n_u_onsite_s_e_' + ":" + repr((i, j, l, n)))
    # bm.jl:2061 — constraint n_u_ehf_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for f in Fu_n[l]:
                        for n in Nd:
                            if n_filter[l, k]  and  nd_en_fltr[n]:
                                m.add((n_u_ehf_d_[i, j, l, k, f, n]) - (( p.n_c_Hr[f, l, k, n]*p.n_c_Fgenf[f, l, k, n] *n_u_onsite_d_[i, j, l, k, n] )), '==', 0, name='n_u_ehf_e_' + ":" + repr((i, j, l, k, f, n)))
    # bm.jl:2071 — constraint n_u_ehf_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for f in Fu_n[l]:
                        for n in Nd:
                            if n_filter[l, k]  and  nd_en_fltr[n]:
                                m.add((n_u_ehf_d_[i, j, l, k, f, n]) - (p.n_u_ehf_ub[f, l, n] * y_n[i, j, l, k]), '<=', 0, name='n_u_ehf_i_' + ":" + repr((i, j, l, k, f, n)))
    # bm.jl:2077 — constraint n_u_ehf_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for f in Fu_n[l]:
                        for n in Nd:
                            if n_filter[l, k]  and  nd_en_fltr[n]:
                                m.add((n_u_ehf[i, j, l, f, n]) - (sum(n_u_ehf_d_[i, j, l, k, f, n] for k in Kn if n_filter[l, k])), '==', 0, name='n_u_ehf_s_e_' + ":" + repr((i, j, l, k, f, n)))
    # bm.jl:2086 — constraint n_cpe_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if n_filter[l, k]  and  nd_em_fltr[n]:
                            m.add((n_cpe_d_[i, j, l, k, n]) - (p.n_c_cpe[l, k, n] * n_cp_d_[i, j, l, k, n] + p.n_rhs_cpe[l, k, n] * y_n[i, j, l, k]), '==', 0, name='n_cpe_d_e_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:2093 — constraint n_cpe_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if n_filter[l, k]  and  nd_em_fltr[n]:
                            m.add((n_cpe_d_[i, j, l, k, n]) - (p.n_cpe_ub[l, n] * y_n[i, j, l, k]), '<=', 0, name='n_cpe_m_i_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:2097 — constraint n_cpe_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_em_fltr[n]:
                        m.add((n_cpe[i, j, l, n]) - (sum(n_cpe_d_[i, j, l, k, n] for k in Kn if n_filter[l, k])), '==', 0, name='n_cpe_s_e_' + ":" + repr((i, j, l, n)))
    # bm.jl:2102 — constraint n_fu_e_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if n_filter[l, k]  and  nd_en_fltr[n]:
                            m.add((n_fu_e_d_[i, j, l, k, n]) - (sum(p.n_c_Fe[f, l, k, n] * n_ehf_d_[i, j, l, k, f, n] for f in Fu_n[l])), '==', 0, name='n_fu_e_d_e_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:2107 — constraint n_fu_e_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if n_filter[l, k]  and  nd_en_fltr[n]:
                            m.add((n_fu_e_d_[i, j, l, k, n]) - (p.n_fu_e_ub[l, n] * y_n[i, j, l, k]), '<=', 0, name='n_fu_e_m_i_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:2112 — constraint n_fu_e_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        m.add((n_fu_e[i, j, l, n]) - (sum(n_fu_e_d_[i, j, l, k, n] for k in Kn if n_filter[l, k])), '==', 0, name='n_fu_e_s_e_' + ":" + repr((i, j, l, n)))
    # bm.jl:2117 — constraint n_u_fu_e_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if n_filter[l, k]  and  nd_en_fltr[n]:
                            m.add((n_u_fu_e_d_[i, j, l, k, n]) - (sum(p.n_c_Fe[f, l, k, n] * n_u_ehf_d_[i, j, l, k, f, n] for f in Fu_n[l])), '==', 0, name='n_u_fu_e_d_e_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:2123 — constraint n_u_fu_e_d_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        if n_filter[l, k]  and  nd_en_fltr[n]:
                            m.add((n_u_fu_e_d_[i, j, l, k, n]) - (p.n_u_fu_e_ub[l, n] * y_n[i, j, l, k]), '<=', 0, name='n_u_fu_e_d_m_i_' + ":" + repr((i, j, l, k, n)))
    # bm.jl:2128 — constraint n_u_fu_e_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        m.add((n_u_fu_e[i, j, l, n]) - (sum(n_u_fu_e_d_[i, j, l, k, n] for k in Kn if n_filter[l, k])), '==', 0, name='n_u_fu_e_s_e_' + ":" + repr((i, j, l, n)))
    # bm.jl:2138 — constraint n_ep0_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        m.add((n_ep0_d_[i, j, l, k]) - (sum(n_fu_e_d_[i, j, l, k, n] for n in Nd if nd_en_fltr[n]) + sum(n_cpe_d_[i, j, l, k, n] for n in Nd if nd_em_fltr[n]) + sum(n_u_fu_e_d_[i, j, l, k, n] for n in Nd if nd_en_fltr[n])), '==', 0, name='n_ep0_d_e_' + ":" + repr((i, j, l, k)))
    # bm.jl:2148 — constraint n_ep0_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        m.add((n_ep0_d_[i, j, l, k]) - (p.n_ep0_bM[l] * y_n[i, j, l, k]), '<=', 0, name='n_ep0_m_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:2151 — constraint n_ep0_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((n_ep0[i, j, l]) - (sum(n_ep0_d_[i, j, l, k] for k in Kn if n_filter[l, k])), '==', 0, name='n_ep0_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:2157 — constraint n_ep1ge_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        m.add((n_ep1ge_d_[i, j, l, k]) - ((1.0 - p.n_chi[l, k])*n_ep0_d_[i, j, l, k]), '==', 0, name='n_ep1ge_d_e_' + ":" + repr((i, j, l, k)))
    # bm.jl:2162 — constraint n_ep1ge_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        m.add((n_ep1ge_d_[i, j, l, k]) - (p.n_ep1ge_bM[l] * y_n[i, j, l, k]), '<=', 0, name='n_ep1ge_m_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:2165 — constraint n_ep1ge_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((n_ep1ge[i, j, l]) - (sum(n_ep1ge_d_[i, j, l, k] for k in Kn if n_filter[l, k])), '==', 0, name='n_ep1ge_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:2172 — constraint n_ep1gce_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        m.add((n_ep1gce_d_[i, j, l, k]) - (p.n_chi[l, k] * (1 - p.n_sigma[l, k]) * n_ep0_d_[i, j, l, k]), '==', 0, name='n_ep1gce_d_e_' + ":" + repr((i, j, l, k)))
    # bm.jl:2177 — constraint n_ep1gce_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        m.add((n_ep1gce_d_[i, j, l, k]) - (p.n_ep1gce_bM[l] * y_n[i, j, l, k]), '<=', 0, name='n_ep1gce_m_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:2180 — constraint n_ep1gce_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((n_ep1gce[i, j, l]) - (sum(n_ep1gce_d_[i, j, l, k] for k in Kn if n_filter[l, k])), '==', 0, name='n_ep1gce_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:2186 — constraint n_ep1gcs_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        m.add((n_ep1gcs_d_[i, j, l, k]) - (p.n_chi[l, k] * p.n_sigma[l, k] * n_ep0_d_[i, j, l, k]), '==', 0, name='n_ep1gcs_d_e_' + ":" + repr((i, j, l, k)))
    # bm.jl:2191 — constraint n_ep1gcs_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        m.add((n_ep1gcs_d_[i, j, l, k]) - (p.n_ep1gcs_bM[l] * y_n[i, j, l, k]), '<=', 0, name='n_ep1gcs_m_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:2194 — constraint n_ep1gcs_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((n_ep1gcs[i, j, l]) - (sum(n_ep1gcs_d_[i, j, l, k] for k in Kn if n_filter[l, k])), '==', 0, name='n_ep1gcs_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:2200 — constraint n_ups_e_mt_in_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        m.add((n_ups_e_mt_in_d_[i, j, l, k]) - (sum( sum(p.n_c_upsein_rate[l, c] * n_x_in_d_[i, j, l, k, n, c] for c in Mt if node_mat[n, c]  and  c  in  input_mat[n] ) for n in Nd)), '==', 0, name='n_ups_e_mt_in_e_' + ":" + repr((i, j, l, k)))
    # bm.jl:2209 — constraint n_ups_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        m.add((n_ups_e_mt_in_d_[i, j, l, k]) - (p.n_ups_e_mt_in_ub[l] * y_n[i, j, l, k]), '<=', 0, name='n_ups_m_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:2212 — constraint n_ups_e_mt_in_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((n_ups_e_mt_in[i, j, l]) - (sum(n_ups_e_mt_in_d_[i, j, l, k] for k in Kn if n_filter[l, k])), '==', 0, name='n_ups_e_mt_in_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:2250 — constraint n_cfonm_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        m.add((n_cfonm_d_[i, j, l, k]) - ((p.n_c_fOnm[l, k] * n_c0_d_[i, j, l, k] + p.n_rhs_fOnm[l, k] * y_n[i, j, l, k])), '==', 0, name='n_cfonm_d_e_' + ":" + repr((i, j, l, k)))
    # bm.jl:2255 — constraint n_cfonm_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        m.add((n_cfonm_d_[i, j, l, k]) - (p.n_cfonm_bM[l] * y_n[i, j, l, k]), '<=', 0, name='n_cfonm_m_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:2258 — constraint n_cfonm_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((n_cfonm[i, j, l]) - (sum(n_cfonm_d_[i, j, l, k] for k in Kn if n_filter[l, k])), '==', 0, name='n_cfonm_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:2264 — constraint n_cvonm_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        m.add((n_cvonm_d_[i, j, l, k]) - ((p.n_c_vOnm[l, k] * n_cp_d_[i, j, l, k, key_node] + p.n_rhs_vOnm[l, k] * y_n[i, j, l, k])), '==', 0, name='n_cvonm_d_e_' + ":" + repr((i, j, l, k)))
    # bm.jl:2269 — constraint n_cvonm_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if n_filter[l, k]:
                        m.add((n_cvonm_d_[i, j, l, k]) - (p.n_cvonm_bM[l] * y_n[i, j, l, k]), '<=', 0, name='n_cvonm_m_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:2272 — constraint n_cvonm_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((n_cvonm[i, j, l]) - (sum(n_cvonm_d_[i, j, l, k] for k in Kn if n_filter[l, k])), '==', 0, name='n_cvonm_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:2281 — constraint n_logic_0_y_i_
    for i in P:
        for j in P2:
            for l in L:
                if j<n_subperiods:
                    m.add((y_n[i, j+1, l, fKn]) - (y_n[i, j, l, fKn]), '<=', 0, name='n_logic_0_y_i_' + ":" + repr((i, j, l)))
    # bm.jl:2284 — constraint n_logic_1_y_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if j<n_subperiods  and  k>fKn  and  n_filter[l, k]:
                        m.add((y_n[i, j+1, l, k]) - (y_n[i, j, l, k]), '>=', 0, name='n_logic_1_y_i_' + ":" + repr((i, j, l, k)))
    # bm.jl:2288 — constraint n_logic_s_e
    for i in P:
        for j in P2:
            for l in L:
                m.add((sum(y_n[i, j, l, k] for k in Kn if n_filter[l, k])) - (1), '==', 0, name='n_logic_s_e' + ":" + repr((i, j, l)))
    # bm.jl:2293 — constraint n_logic_2
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    if k>fKn  and  (i, j)!=(fP,fP2)  and  n_filter[l, k]:
                        m.add((y_n[i, j, l, k] + y_o[i, j, l]) - (1), '<=', 0, name='n_logic_2' + ":" + repr((i, j, l, k)))
    # bm.jl:2302 — constraint n_u_cost_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((n_u_cost[i, j, l]) - (sum(p.c_u_cost[i, j, l] * n_u[i, j, l, n] for n in Nd if nd_en_fltr[n])), '==', 0, name='n_u_cost_e_' + ":" + repr((i, j, l)))
    # bm.jl:2306 — constraint n_ehf_cost_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((n_ehf_cost[i, j, l]) - (sum(p.c_n_ehf_cost[i, j, l, f]* (n_ehf[i, j, l, f, n] + n_u_ehf[i, j, l, f, n]) for f in Fu_n[l] for n in Nd if nd_en_fltr[n])), '==', 0, name='n_ehf_cost_e_' + ":" + repr((i, j, l)))
    # bm.jl:2315 — constraint o_cp_d_e
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    m.add((o_cp_d_[i, j, l, n]) - (o_tcp_d_[i, j, l, sT, n]), '==', 0, name='o_cp_d_e' + ":" + repr((i, j, l, n)))
    # bm.jl:2318 — constraint o_cp_m1_i_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    m.add((o_cp_d_[i, j, l, n]) - (p.o_cp_ub[l, n] * y_o[i, j, l]), '<=', 0, name='o_cp_m1_i_' + ":" + repr((i, j, l, n)))
    # bm.jl:2322 — constraint o_cp_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    m.add((o_cp[i, j, l, n]) - (o_cp_d_[i, j, l, n]), '==', 0, name='o_cp_s_e_' + ":" + repr((i, j, l, n)))
    # bm.jl:2326 — constraint o_tcp_d_m1_i_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    m.add((o_tcp_d_[i, j, l, sT, n]) - (p.o_cp_ub[l, n] * y_o[i, j, l]), '<=', 0, name='o_tcp_d_m1_i_' + ":" + repr((i, j, l, n)))
    # bm.jl:2329 — constraint o_tcp_d_m0_i_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    m.add((o_tcp_d_[i, j, l, sF, n]) - (p.o_cp_ub[l, n] * (1 - y_o[i, j, l])), '<=', 0, name='o_tcp_d_m0_i_' + ":" + repr((i, j, l, n)))
    # bm.jl:2332 — constraint o_tcp_d_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    m.add((r_cp[i, j, l, n]) - (o_tcp_d_[i, j, l, sT, n] + o_tcp_d_[i, j, l, sF, n]), '==', 0, name='o_tcp_d_s_e_' + ":" + repr((i, j, l, n)))
    # bm.jl:2357 — constraint n_mass_b_e_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        for c in Mt:
                            if n_filter[l, k]  and  node_mat[n, c]  and  not skip_mb[n, c]:
                                m.add((f_xin(m, i, j, l, k, n, c, input_mat, mode="n") - f_xout(m, i, j, l, k, n, c, output_mat, mode="n") + (p.n_Kmb[n, c, ckey[n], k] * n_cp_d_[i, j, l, k, n])) - (0.0), '==', 0, name='n_mass_b_e_' + ":" + repr((i, j, l, k, n, c)))
    # bm.jl:2366 — constraint n_comp_link_e_
    for i in P:
        for j in P2:
            for l in L:
                for lnk in Ln:
                    m.add((n_x_out[i, j, l, links_list[lnk][1], links_list[lnk][3]]) - (n_x_in[i, j, l, links_list[lnk][2], links_list[lnk][3]]), '==', 0, name='n_comp_link_e_' + ":" + repr((i, j, l, lnk)))
    # bm.jl:2372 — constraint n_x_in_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        for c in Mt:
                            if n_filter[l, k]  and  node_mat[n, c]  and  c  in  input_mat[n]:
                                m.add((n_x_in_d_[i, j, l, k, n, c]) - (p.n_x_in_ub[l, n, c] * y_n[i, j, l, k]), '<=', 0, name='n_x_in_m_i_' + ":" + repr((i, j, l, k, n, c)))
    # bm.jl:2378 — constraint n_x_in_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  input_mat[n]:
                            m.add((n_x_in[i, j, l, n, c]) - (sum(n_x_in_d_[i, j, l, k, n, c] for k in Kn if n_filter[l, k])), '==', 0, name='n_x_in_s_e_' + ":" + repr((i, j, l, n, c)))
    # bm.jl:2384 — constraint n_x_out_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for k in Kn:
                    for n in Nd:
                        for c in Mt:
                            if n_filter[l, k]  and  node_mat[n, c]  and  c  in  output_mat[n]:
                                m.add((n_x_out_d_[i, j, l, k, n, c]) - (p.n_x_out_ub[l, n, c] * y_n[i, j, l, k]), '<=', 0, name='n_x_out_m_i_' + ":" + repr((i, j, l, k, n, c)))
    # bm.jl:2390 — constraint n_x_out_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    for c in Mt:
                        if node_mat[n, c]  and  c  in  output_mat[n]:
                            m.add((n_x_out[i, j, l, n, c]) - (sum(n_x_out_d_[i, j, l, k, n, c] for k in Kn if n_filter[l, k])), '==', 0, name='n_x_out_s_e_' + ":" + repr((i, j, l, n, c)))
    # bm.jl:2400 — constraint o_u_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        m.add((o_u_d_[i, j, l, n]) - (o_tu_d_[i, j, l, sT, n]), '==', 0, name='o_u_d_e_' + ":" + repr((i, j, l, n)))
    # bm.jl:2404 — constraint o_u_m_i_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        m.add((o_u_d_[i, j, l, n]) - (p.o_u_ub[l, n] * y_o[i, j, l]), '<=', 0, name='o_u_m_i_' + ":" + repr((i, j, l, n)))
    # bm.jl:2408 — constraint o_u_s_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        m.add((o_u[i, j, l, n]) - (o_u_d_[i, j, l, n]), '==', 0, name='o_u_s_' + ":" + repr((i, j, l, n)))
    # bm.jl:2412 — constraint o_tu_d_m1_i_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        m.add((o_tu_d_[i, j, l, sT, n]) - (p.o_u_ub[l, n] * y_o[i, j, l]), '<=', 0, name='o_tu_d_m1_i_' + ":" + repr((i, j, l, n)))
    # bm.jl:2416 — constraint o_tu_d_m0_i_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        m.add((o_tu_d_[i, j, l, sF, n]) - (p.o_u_ub[l, n] * (1 - y_o[i, j, l])), '<=', 0, name='o_tu_d_m0_i_' + ":" + repr((i, j, l, n)))
    # bm.jl:2420 — constraint o_tu_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for n in Nd:
                    if nd_en_fltr[n]:
                        m.add((r_u[i, j, l, n]) - (o_tu_d_[i, j, l, sT, n] + o_tu_d_[i, j, l, sF, n]), '==', 0, name='o_tu_s_e_' + ":" + repr((i, j, l, n)))
    # bm.jl:2429 — constraint o_ehf_d_e_
    for i in P:
        for j in P2:
            for l in L:
                for f in Fu_r[l]:
                    for n in Nd:
                        if nd_en_fltr[n]:
                            m.add((o_ehf_d_[i, j, l, f, n]) - (o_tehf_d_[i, j, l, f, sT, n]), '==', 0, name='o_ehf_d_e_' + ":" + repr((i, j, l, f, n)))
    # bm.jl:2433 — constraint o_ehf_m_i
    for i in P:
        for j in P2:
            for l in L:
                for f in Fu_r[l]:
                    for n in Nd:
                        if nd_en_fltr[n]:
                            m.add((o_ehf_d_[i, j, l, f, n]) - (p.o_ehf_ub[f, l, n] * y_o[i, j, l]), '<=', 0, name='o_ehf_m_i' + ":" + repr((i, j, l, f, n)))
    # bm.jl:2437 — constraint o_ehf_s_
    for i in P:
        for j in P2:
            for l in L:
                for f in Fu_r[l]:
                    for n in Nd:
                        if nd_en_fltr[n]:
                            m.add((o_ehf[i, j, l, f, n]) - (o_ehf_d_[i, j, l, f, n]), '==', 0, name='o_ehf_s_' + ":" + repr((i, j, l, f, n)))
    # bm.jl:2441 — constraint o_tehf_d_m1_i_
    for i in P:
        for j in P2:
            for l in L:
                for f in Fu_r[l]:
                    for n in Nd:
                        if nd_en_fltr[n]:
                            m.add((o_tehf_d_[i, j, l, f, sT, n]) - (p.o_ehf_ub[f, l, n] * y_o[i, j, l]), '<=', 0, name='o_tehf_d_m1_i_' + ":" + repr((i, j, l, f, n)))
    # bm.jl:2445 — constraint o_tehf_d_m0_i_
    for i in P:
        for j in P2:
            for l in L:
                for f in Fu_r[l]:
                    for n in Nd:
                        if nd_en_fltr[n]:
                            m.add((o_tehf_d_[i, j, l, f, sF, n]) - (p.o_ehf_ub[f, l, n] * (1 - y_o[i, j, l])), '<=', 0, name='o_tehf_d_m0_i_' + ":" + repr((i, j, l, f, n)))
    # bm.jl:2450 — constraint o_tehf_s_e_
    for i in P:
        for j in P2:
            for l in L:
                for f in Fu_r[l]:
                    for n in Nd:
                        if nd_en_fltr[n]:
                            m.add((r_ehf[i, j, l, f, n] + r_u_ehf[i, j, l, f, n]) - (o_tehf_d_[i, j, l, f, sT, n] + o_tehf_d_[i, j, l, f, sF, n]), '==', 0, name='o_tehf_s_e_' + ":" + repr((i, j, l, f, n)))
    # bm.jl:2458 — constraint o_ep0_d_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_ep0_d_[i, j, l]) - (o_tep0_d_[i, j, l, sT]), '==', 0, name='o_ep0_d_e_' + ":" + repr((i, j, l)))
    # bm.jl:2461 — constraint o_ep0_m_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_ep0_d_[i, j, l]) - (p.o_ep0_bM[l] * y_o[i, j, l]), '<=', 0, name='o_ep0_m_i_' + ":" + repr((i, j, l)))
    # bm.jl:2464 — constraint o_ep0_s_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_ep0[i, j, l]) - (o_ep0_d_[i, j, l]), '==', 0, name='o_ep0_s_' + ":" + repr((i, j, l)))
    # bm.jl:2467 — constraint o_tep0_d_m1_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_tep0_d_[i, j, l, sT]) - (p.o_ep0_bM[l] * y_o[i, j, l]), '<=', 0, name='o_tep0_d_m1_i_' + ":" + repr((i, j, l)))
    # bm.jl:2470 — constraint o_tep0_d_m0_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_tep0_d_[i, j, l, sF]) - (p.o_ep0_bM[l] * (1 - y_o[i, j, l])), '<=', 0, name='o_tep0_d_m0_i_' + ":" + repr((i, j, l)))
    # bm.jl:2473 — constraint o_tep0_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_ep0[i, j, l]) - (o_tep0_d_[i, j, l, sT] + o_tep0_d_[i, j, l, sF]), '==', 0, name='o_tep0_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:2478 — constraint o_ep1ge_d_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_ep1ge_d_[i, j, l]) - (o_tep1ge_d_[i, j, l, sT]), '==', 0, name='o_ep1ge_d_e_' + ":" + repr((i, j, l)))
    # bm.jl:2481 — constraint o_ep1ge_m_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_ep1ge_d_[i, j, l]) - (p.o_ep1ge_bM[l] * y_o[i, j, l]), '<=', 0, name='o_ep1ge_m_i_' + ":" + repr((i, j, l)))
    # bm.jl:2484 — constraint o_ep1ge_s_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_ep1ge[i, j, l]) - (o_ep1ge_d_[i, j, l]), '==', 0, name='o_ep1ge_s_' + ":" + repr((i, j, l)))
    # bm.jl:2487 — constraint o_tep1ge_d_m1_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_tep1ge_d_[i, j, l, sT]) - (p.o_ep1ge_bM[l] * y_o[i, j, l]), '<=', 0, name='o_tep1ge_d_m1_i_' + ":" + repr((i, j, l)))
    # bm.jl:2490 — constraint o_tep1ge_d_m0_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_tep1ge_d_[i, j, l, sF]) - (p.o_ep1ge_bM[l] * (1 - y_o[i, j, l])), '<=', 0, name='o_tep1ge_d_m0_i_' + ":" + repr((i, j, l)))
    # bm.jl:2493 — constraint o_tep1ge_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_ep1ge[i, j, l]) - (o_tep1ge_d_[i, j, l, sT] + o_tep1ge_d_[i, j, l, sF]), '==', 0, name='o_tep1ge_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:2499 — constraint o_ep1gce_d_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_ep1gce_d_[i, j, l]) - (o_tep1gce_d_[i, j, l, sT]), '==', 0, name='o_ep1gce_d_e_' + ":" + repr((i, j, l)))
    # bm.jl:2502 — constraint o_ep1gce_m_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_ep1gce_d_[i, j, l]) - (p.o_ep1gce_bM[l] * y_o[i, j, l]), '<=', 0, name='o_ep1gce_m_i_' + ":" + repr((i, j, l)))
    # bm.jl:2505 — constraint o_ep1gce_s_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_ep1gce[i, j, l]) - (o_ep1gce_d_[i, j, l]), '==', 0, name='o_ep1gce_s_' + ":" + repr((i, j, l)))
    # bm.jl:2508 — constraint o_tep1gce_d_m1_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_tep1gce_d_[i, j, l, sT]) - (p.o_ep1gce_bM[l] * y_o[i, j, l]), '<=', 0, name='o_tep1gce_d_m1_i_' + ":" + repr((i, j, l)))
    # bm.jl:2511 — constraint o_tep1gce_d_m0_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_tep1gce_d_[i, j, l, sT]) - (p.o_ep1gce_bM[l] * (1 - y_o[i, j, l])), '<=', 0, name='o_tep1gce_d_m0_i_' + ":" + repr((i, j, l)))
    # bm.jl:2514 — constraint o_tep1gce_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_ep1gce[i, j, l]) - (o_tep1gce_d_[i, j, l, sT] + o_tep1gce_d_[i, j, l, sF]), '==', 0, name='o_tep1gce_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:2521 — constraint o_ep1gcs_d_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_ep1gcs_d_[i, j, l]) - (o_tep1gcs_d_[i, j, l, sT]), '==', 0, name='o_ep1gcs_d_e_' + ":" + repr((i, j, l)))
    # bm.jl:2524 — constraint o_ep1gcs_m_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_ep1gcs_d_[i, j, l]) - (p.o_ep1gcs_bM[l] * y_o[i, j, l]), '<=', 0, name='o_ep1gcs_m_i_' + ":" + repr((i, j, l)))
    # bm.jl:2527 — constraint o_ep1gcs_s_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_ep1gcs[i, j, l]) - (o_ep1gcs_d_[i, j, l]), '==', 0, name='o_ep1gcs_s_' + ":" + repr((i, j, l)))
    # bm.jl:2530 — constraint o_tep1gcs_d_m1_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_tep1gcs_d_[i, j, l, sT]) - (p.o_ep1gcs_bM[l] * y_o[i, j, l]), '<=', 0, name='o_tep1gcs_d_m1_i_' + ":" + repr((i, j, l)))
    # bm.jl:2533 — constraint o_tepgcs_d_m0_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_tep1gcs_d_[i, j, l, sF]) - (p.o_ep1gcs_bM[l] * (1 - y_o[i, j, l])), '<=', 0, name='o_tepgcs_d_m0_i_' + ":" + repr((i, j, l)))
    # bm.jl:2536 — constraint o_tep1gcs_s_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((r_ep1gcs[i, j, l]) - (o_tep1gcs_d_[i, j, l, sT] + o_tep1gcs_d_[i, j, l, sF]), '==', 0, name='o_tep1gcs_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:2572 — constraint o_last_loan_s_e_
    for i in P:
        for l in L:
            m.add((sum(r_loan[i, n_subperiods, i0, j0, l] for j0 in P2 for i0 in P if _in((i, n_subperiods),loan_window[i0, j0])  and  _ne((i,n_subperiods),(i0,j0)) ) + sum(e_loan[i, n_subperiods, i0, j0, l] for j0 in P2 for i0 in P if _in((i, n_subperiods),loan_window[i0, j0])  and  _ne((i,n_subperiods),(i0,j0)) )) - (o_loan_last[i, l, sT] + o_loan_last[i, l, sF]), '==', 0, name='o_last_loan_s_e_' + ":" + repr((i, l)))
    # bm.jl:2582 — constraint o_last_loan_d_m1_i_
    for i in P:
        for l in L:
            m.add((o_loan_last[i, l, sT]) - (oloanlastub[l] * y_o[i, n_subperiods, l]), '<=', 0, name='o_last_loan_d_m1_i_' + ":" + repr((i, l)))
    # bm.jl:2585 — constraint o_last_loan_d_m0_i_
    for i in P:
        for l in L:
            m.add((o_loan_last[i, l, sF]) - (oloanlastub[l] * (1 - y_o[i, n_subperiods, l])), '<=', 0, name='o_last_loan_d_m0_i_' + ":" + repr((i, l)))
    # bm.jl:2590 — constraint o_u_cost_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_u_cost[i, j, l]) - (sum( p.c_u_cost[i, j, l] * o_u[i, j, l, n] for n in Nd if nd_en_fltr[n] )), '==', 0, name='o_u_cost_e_' + ":" + repr((i, j, l)))
    # bm.jl:2597 — constraint o_ehf_cost_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_ehf_cost[i, j, l]) - (sum( sum(p.c_r_ehf_cost[i, j, l, f] * o_ehf[i, j, l, f, n] for f in Fu_r[l]) for n in Nd if nd_en_fltr[n] )), '==', 0, name='o_ehf_cost_e_' + ":" + repr((i, j, l)))
    # bm.jl:2605 — constraint o_ep1gcs_cost_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_ep1gcs_cost[i, j, l]) - (p.c_cts_cost[l] * o_ep1gcs[i, j, l]), '==', 0, name='o_ep1gcs_cost_e_' + ":" + repr((i, j, l)))
    # bm.jl:2609 — constraint n_ep1gcs_cost_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((n_ep1gcs_cost[i, j, l]) - (p.c_cts_cost[l] * n_ep1gcs[i, j, l]), '==', 0, name='n_ep1gcs_cost_e_' + ":" + repr((i, j, l)))
    # bm.jl:2614 — constraint o_x_in_cost_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_x_in_cost[i, j, l]) - (sum(p.c_xin_cost[l, c] * o_x_in[i, j, l, n, c] for n in Nd for c in Mt if node_mat[n, c]  and  c  in  input_mat[n] )), '==', 0, name='o_x_in_cost_e_' + ":" + repr((i, j, l)))
    # bm.jl:2620 — constraint n_x_in_cost_e_
    for i in P:
        for j in P2:
            for l in L:
                m.add((n_x_in_cost[i, j, l]) - (sum(p.c_xin_cost[l, c] * n_x_in[i, j, l, n, c] for n in Nd for c in Mt if node_mat[n, c]  and  c  in  input_mat[n] )), '==', 0, name='n_x_in_cost_e_' + ":" + repr((i, j, l)))
    attachInitCond(index_l, m, p, s) if not isinstance(index_p, int) or index_p == fP else None
    return m

def attachInitCond(index_l, m, p, s):
    P = s.P
    P2 = s.P2
    L = s.L
    Kr = s.Kr
    Kn = s.Kn
    Fu_r = s.Fu_r
    Fu_n = s.Fu_n
    Nf = s.Nf
    Nd = s.Nd
    Mt = s.Mt
    Ln = s.Ln
    n_periods = p.n_periods
    n_subperiods = p.n_subperiods
    yr_subperiod = p.yr_subperiod
    loan_subperiods = 30 // yr_subperiod
    loan_periods = loan_subperiods // n_subperiods
    remain_loan_years = 30 - loan_subperiods * yr_subperiod
    remain_subp_loan = loan_subperiods - loan_periods * n_subperiods
    remain_sub_loan = remain_subp_loan + int(remain_loan_years > 0)
    sT = 1
    sF = 2
    fP = min(P)
    fP2 = min(P2)
    fKr = min(Kr)
    fKn = min(Kn)
    r_filter = p.r_filter
    n_filter = p.n_filter
    nd_en_fltr = p.nd_en_fltr
    nd_em_fltr = p.nd_em_fltr
    node_mat = p.node_mat
    skip_mb = p.skip_mb
    links_list = p.links_list
    input_mat = p.input_mat
    output_mat = p.output_mat
    ckey = p.ckey
    key_node = p.key_node
    loan_window, pay_window = gen_lp_window(p, s)
    y_o = m.y_o
    y_r = m.y_r
    y_e = m.y_e
    e_c_d_ = m.e_c_d_
    e_c = m.e_c
    e_l = m.e_l
    x = m.x
    x_d_ = m.x_d_
    r_cp = m.r_cp
    r_cp_d_ = m.r_cp_d_
    cpb = m.cpb
    r_cpb_d_ = m.r_cpb_d_
    r_eh = m.r_eh
    r_eh_d_ = m.r_eh_d_
    r_ehf = m.r_ehf
    r_ehf_d_ = m.r_ehf_d_
    r_u = m.r_u
    r_u_d_ = m.r_u_d_
    r_u_onsite = m.r_u_onsite
    r_u_onsite_d_ = m.r_u_onsite_d_
    r_u_ehf = m.r_u_ehf
    r_u_ehf_d_ = m.r_u_ehf_d_
    r_cpe = m.r_cpe
    r_cpe_d_ = m.r_cpe_d_
    r_fu_e = m.r_fu_e
    r_fu_e_d_ = m.r_fu_e_d_
    r_u_fu_e = m.r_u_fu_e
    r_u_fu_e_d_ = m.r_u_fu_e_d_
    r_ep0 = m.r_ep0
    r_ep0_d_ = m.r_ep0_d_
    r_ep1ge = m.r_ep1ge
    r_ep1ge_d_ = m.r_ep1ge_d_
    r_ep1gce = m.r_ep1gce
    r_ep1gce_d_ = m.r_ep1gce_d_
    r_ep1gcs = m.r_ep1gcs
    r_ep1gcs_d_ = m.r_ep1gcs_d_
    r_fstck = m.r_fstck
    r_fstck_d_ = m.r_fstck_d_
    r_cfonm_d_ = m.r_cfonm_d_
    r_cfonm = m.r_cfonm
    r_cvonm_d_ = m.r_cvonm_d_
    r_cvonm = m.r_cvonm
    r_loan = m.r_loan
    r_l0 = m.r_l0
    r_l0_d_ = m.r_l0_d_
    r_l0_pd_ = m.r_l0_pd_
    r_l0add = m.r_l0add
    r_l0add_d_ = m.r_l0add_d_
    r_e_c_d_ = m.r_e_c_d_
    r_le = m.r_le
    r_le_d_ = m.r_le_d_
    r_le_pd_ = m.r_le_pd_
    r_leadd = m.r_leadd
    r_leadd_d_ = m.r_leadd_d_
    r_le_ped_ = m.r_le_ped_
    r_leadde = m.r_leadde
    r_leadde_d_ = m.r_leadde_d_
    r_pay = m.r_pay
    r_paye = m.r_paye
    e_ladd = m.e_ladd
    e_ladd_d_ = m.e_ladd_d_
    e_l_d_ = m.e_l_d_
    e_l_pd_ = m.e_l_pd_
    e_ann = m.e_ann
    e_ann_d_ = m.e_ann_d_
    e_ann_0 = m.e_ann_0
    e_ann_1 = m.e_ann_1
    e_pay = m.e_pay
    e_loan = m.e_loan
    t_loan_d_ = m.t_loan_d_
    t_ret_cost = m.t_ret_cost
    t_ret_cost_d_ = m.t_ret_cost_d_
    y_n = m.y_n
    n_c0 = m.n_c0
    n_c0_d_ = m.n_c0_d_
    n_cp = m.n_cp
    n_cp_d_ = m.n_cp_d_
    n_l_d_ = m.n_l_d_
    n_l = m.n_l
    n_ann_d_ = m.n_ann_d_
    n_ann = m.n_ann
    n_ann_0 = m.n_ann_0
    n_ann_1 = m.n_ann_1
    n_ladd_d_ = m.n_ladd_d_
    n_ladd = m.n_ladd
    n_l_pd_ = m.n_l_pd_
    n_loan = m.n_loan
    n_pay = m.n_pay
    n_eh = m.n_eh
    n_eh_d_ = m.n_eh_d_
    n_ehf = m.n_ehf
    n_ehf_d_ = m.n_ehf_d_
    n_u = m.n_u
    n_u_d_ = m.n_u_d_
    n_u_onsite = m.n_u_onsite
    n_u_onsite_d_ = m.n_u_onsite_d_
    n_u_ehf = m.n_u_ehf
    n_u_ehf_d_ = m.n_u_ehf_d_
    n_cpe = m.n_cpe
    n_cpe_d_ = m.n_cpe_d_
    n_fu_e = m.n_fu_e
    n_fu_e_d_ = m.n_fu_e_d_
    n_u_fu_e = m.n_u_fu_e
    n_u_fu_e_d_ = m.n_u_fu_e_d_
    n_ep0 = m.n_ep0
    n_ep0_d_ = m.n_ep0_d_
    n_ep1ge = m.n_ep1ge
    n_ep1ge_d_ = m.n_ep1ge_d_
    n_ep1gce = m.n_ep1gce
    n_ep1gce_d_ = m.n_ep1gce_d_
    n_ep1gcs = m.n_ep1gcs
    n_ep1gcs_d_ = m.n_ep1gcs_d_
    n_fstck = m.n_fstck
    n_fstck_d_ = m.n_fstck_d_
    n_cfonm_d_ = m.n_cfonm_d_
    n_cfonm = m.n_cfonm
    n_cvonm_d_ = m.n_cvonm_d_
    n_cvonm = m.n_cvonm
    n_u_cost = m.n_u_cost
    n_ehf_cost = m.n_ehf_cost
    o_pay = m.o_pay
    o_pay_d_ = m.o_pay_d_
    o_cfonm = m.o_cfonm
    o_cfonm_d_ = m.o_cfonm_d_
    o_tcfonm_d_ = m.o_tcfonm_d_
    o_cvonm = m.o_cvonm
    o_cvonm_d_ = m.o_cvonm_d_
    o_tcvonm_d_ = m.o_tcvonm_d_
    o_cp = m.o_cp
    o_cp_d_ = m.o_cp_d_
    o_tcp_d_ = m.o_tcp_d_
    o_tpay_d_ = m.o_tpay_d_
    o_rcp_d_ = m.o_rcp_d_
    o_rcp = m.o_rcp
    o_u_d_ = m.o_u_d_
    o_u = m.o_u
    o_tu_d_ = m.o_tu_d_
    o_ehf_d_ = m.o_ehf_d_
    o_ehf = m.o_ehf
    o_tehf_d_ = m.o_tehf_d_
    o_ep0_d_ = m.o_ep0_d_
    o_ep0 = m.o_ep0
    o_tep0_d_ = m.o_tep0_d_
    o_ep1ge_d_ = m.o_ep1ge_d_
    o_ep1ge = m.o_ep1ge
    o_tep1ge_d_ = m.o_tep1ge_d_
    o_ep1gce_d_ = m.o_ep1gce_d_
    o_ep1gce = m.o_ep1gce
    o_tep1gce_d_ = m.o_tep1gce_d_
    o_ep1gcs_d_ = m.o_ep1gcs_d_
    o_ep1gcs = m.o_ep1gcs
    o_tep1gcs_d_ = m.o_tep1gcs_d_
    o_loan_last = m.o_loan_last
    o_u_cost = m.o_u_cost
    o_ehf_cost = m.o_ehf_cost
    o_ep1gcs_cost = m.o_ep1gcs_cost
    n_ep1gcs_cost = m.n_ep1gcs_cost
    r_x_in_d_ = m.r_x_in_d_
    r_x_out_d_ = m.r_x_out_d_
    r_x_o_s_d_ = m.r_x_o_s_d_
    r_x_in = m.r_x_in
    r_x_out = m.r_x_out
    o_x_in = m.o_x_in
    o_x_in_d_ = m.o_x_in_d_
    o_tx_in_d_ = m.o_tx_in_d_
    o_x_out = m.o_x_out
    o_x_out_d_ = m.o_x_out_d_
    o_tx_out_d_ = m.o_tx_out_d_
    n_x_in_d_ = m.n_x_in_d_
    n_x_out_d_ = m.n_x_out_d_
    n_x_o_s_d_ = m.n_x_o_s_d_
    n_x_in = m.n_x_in
    n_x_out = m.n_x_out
    r_ups_e_mt_in_d_ = m.r_ups_e_mt_in_d_
    r_ups_e_mt_in = m.r_ups_e_mt_in
    o_ups_e_mt_in_d_ = m.o_ups_e_mt_in_d_
    o_ups_e_mt_in = m.o_ups_e_mt_in
    o_tups_e_mt_in_d_ = m.o_tups_e_mt_in_d_
    n_ups_e_mt_in_d_ = m.n_ups_e_mt_in_d_
    n_ups_e_mt_in = m.n_ups_e_mt_in
    o_x_in_cost = m.o_x_in_cost
    n_x_in_cost = m.n_x_in_cost
    L = [index_l] if isinstance(index_l, int) else index_l
    # bm.jl:2688 — constraint r_logic_init_0
    for l in L:
        m.add((y_r[fP, fP2, l, fKr]) - (1), '==', 0, name='r_logic_init_0' + ":" + repr((l,)))
    # bm.jl:2692 — constraint r_logic_init_1
    for l in L:
        for k in Kr:
            if k > fKr  and  r_filter[l, k]:
                m.add((y_r[fP, fP2, l, k]) - (0), '==', 0, name='r_logic_init_1' + ":" + repr((l, k)))
    # bm.jl:2696 — constraint o_logic_init_0_
    for l in L:
        m.add((y_o[fP, fP2, l]) - (1), '==', 0, name='o_logic_init_0_' + ":" + repr((l,)))
    # bm.jl:2699 — constraint n_logic_0_init
    for l in L:
        m.add((y_n[fP, fP2, l, fKn]) - (1), '==', 0, name='n_logic_0_init' + ":" + repr((l,)))
    # bm.jl:2702 — constraint n_logic_1_init
    for l in L:
        for k in Kn:
            if k>fKn  and  n_filter[l, k]:
                m.add((y_n[fP, fP2, l, k]) - (0), '==', 0, name='n_logic_1_init' + ":" + repr((l, k)))
    # bm.jl:2708 — constraint e_ladd_m_0_y_0_
    for l in L:
        m.add((e_ladd_d_[fP, fP2, l]) - (p.e_ladd_ub[l] * y_e[fP, fP2, l]), '<=', 0, name='e_ladd_m_0_y_0_' + ":" + repr((l,)))
    return m

def attachFullObjectiveBlock(m, p, s, tCoeff=1e-6):
    P = s.P
    P2 = s.P2
    L = s.L
    Kr = s.Kr
    Kn = s.Kn
    Fu_r = s.Fu_r
    Fu_n = s.Fu_n
    Nf = s.Nf
    Nd = s.Nd
    Mt = s.Mt
    Ln = s.Ln
    n_periods = p.n_periods
    n_subperiods = p.n_subperiods
    yr_subperiod = p.yr_subperiod
    loan_subperiods = 30 // yr_subperiod
    loan_periods = loan_subperiods // n_subperiods
    remain_loan_years = 30 - loan_subperiods * yr_subperiod
    remain_subp_loan = loan_subperiods - loan_periods * n_subperiods
    remain_sub_loan = remain_subp_loan + int(remain_loan_years > 0)
    sT = 1
    sF = 2
    fP = min(P)
    fP2 = min(P2)
    fKr = min(Kr)
    fKn = min(Kn)
    r_filter = p.r_filter
    n_filter = p.n_filter
    nd_en_fltr = p.nd_en_fltr
    nd_em_fltr = p.nd_em_fltr
    node_mat = p.node_mat
    skip_mb = p.skip_mb
    links_list = p.links_list
    input_mat = p.input_mat
    output_mat = p.output_mat
    ckey = p.ckey
    key_node = p.key_node
    loan_window, pay_window = gen_lp_window(p, s)
    y_o = m.y_o
    y_r = m.y_r
    y_e = m.y_e
    e_c_d_ = m.e_c_d_
    e_c = m.e_c
    e_l = m.e_l
    x = m.x
    x_d_ = m.x_d_
    r_cp = m.r_cp
    r_cp_d_ = m.r_cp_d_
    cpb = m.cpb
    r_cpb_d_ = m.r_cpb_d_
    r_eh = m.r_eh
    r_eh_d_ = m.r_eh_d_
    r_ehf = m.r_ehf
    r_ehf_d_ = m.r_ehf_d_
    r_u = m.r_u
    r_u_d_ = m.r_u_d_
    r_u_onsite = m.r_u_onsite
    r_u_onsite_d_ = m.r_u_onsite_d_
    r_u_ehf = m.r_u_ehf
    r_u_ehf_d_ = m.r_u_ehf_d_
    r_cpe = m.r_cpe
    r_cpe_d_ = m.r_cpe_d_
    r_fu_e = m.r_fu_e
    r_fu_e_d_ = m.r_fu_e_d_
    r_u_fu_e = m.r_u_fu_e
    r_u_fu_e_d_ = m.r_u_fu_e_d_
    r_ep0 = m.r_ep0
    r_ep0_d_ = m.r_ep0_d_
    r_ep1ge = m.r_ep1ge
    r_ep1ge_d_ = m.r_ep1ge_d_
    r_ep1gce = m.r_ep1gce
    r_ep1gce_d_ = m.r_ep1gce_d_
    r_ep1gcs = m.r_ep1gcs
    r_ep1gcs_d_ = m.r_ep1gcs_d_
    r_fstck = m.r_fstck
    r_fstck_d_ = m.r_fstck_d_
    r_cfonm_d_ = m.r_cfonm_d_
    r_cfonm = m.r_cfonm
    r_cvonm_d_ = m.r_cvonm_d_
    r_cvonm = m.r_cvonm
    r_loan = m.r_loan
    r_l0 = m.r_l0
    r_l0_d_ = m.r_l0_d_
    r_l0_pd_ = m.r_l0_pd_
    r_l0add = m.r_l0add
    r_l0add_d_ = m.r_l0add_d_
    r_e_c_d_ = m.r_e_c_d_
    r_le = m.r_le
    r_le_d_ = m.r_le_d_
    r_le_pd_ = m.r_le_pd_
    r_leadd = m.r_leadd
    r_leadd_d_ = m.r_leadd_d_
    r_le_ped_ = m.r_le_ped_
    r_leadde = m.r_leadde
    r_leadde_d_ = m.r_leadde_d_
    r_pay = m.r_pay
    r_paye = m.r_paye
    e_ladd = m.e_ladd
    e_ladd_d_ = m.e_ladd_d_
    e_l_d_ = m.e_l_d_
    e_l_pd_ = m.e_l_pd_
    e_ann = m.e_ann
    e_ann_d_ = m.e_ann_d_
    e_ann_0 = m.e_ann_0
    e_ann_1 = m.e_ann_1
    e_pay = m.e_pay
    e_loan = m.e_loan
    t_loan_d_ = m.t_loan_d_
    t_ret_cost = m.t_ret_cost
    t_ret_cost_d_ = m.t_ret_cost_d_
    y_n = m.y_n
    n_c0 = m.n_c0
    n_c0_d_ = m.n_c0_d_
    n_cp = m.n_cp
    n_cp_d_ = m.n_cp_d_
    n_l_d_ = m.n_l_d_
    n_l = m.n_l
    n_ann_d_ = m.n_ann_d_
    n_ann = m.n_ann
    n_ann_0 = m.n_ann_0
    n_ann_1 = m.n_ann_1
    n_ladd_d_ = m.n_ladd_d_
    n_ladd = m.n_ladd
    n_l_pd_ = m.n_l_pd_
    n_loan = m.n_loan
    n_pay = m.n_pay
    n_eh = m.n_eh
    n_eh_d_ = m.n_eh_d_
    n_ehf = m.n_ehf
    n_ehf_d_ = m.n_ehf_d_
    n_u = m.n_u
    n_u_d_ = m.n_u_d_
    n_u_onsite = m.n_u_onsite
    n_u_onsite_d_ = m.n_u_onsite_d_
    n_u_ehf = m.n_u_ehf
    n_u_ehf_d_ = m.n_u_ehf_d_
    n_cpe = m.n_cpe
    n_cpe_d_ = m.n_cpe_d_
    n_fu_e = m.n_fu_e
    n_fu_e_d_ = m.n_fu_e_d_
    n_u_fu_e = m.n_u_fu_e
    n_u_fu_e_d_ = m.n_u_fu_e_d_
    n_ep0 = m.n_ep0
    n_ep0_d_ = m.n_ep0_d_
    n_ep1ge = m.n_ep1ge
    n_ep1ge_d_ = m.n_ep1ge_d_
    n_ep1gce = m.n_ep1gce
    n_ep1gce_d_ = m.n_ep1gce_d_
    n_ep1gcs = m.n_ep1gcs
    n_ep1gcs_d_ = m.n_ep1gcs_d_
    n_fstck = m.n_fstck
    n_fstck_d_ = m.n_fstck_d_
    n_cfonm_d_ = m.n_cfonm_d_
    n_cfonm = m.n_cfonm
    n_cvonm_d_ = m.n_cvonm_d_
    n_cvonm = m.n_cvonm
    n_u_cost = m.n_u_cost
    n_ehf_cost = m.n_ehf_cost
    o_pay = m.o_pay
    o_pay_d_ = m.o_pay_d_
    o_cfonm = m.o_cfonm
    o_cfonm_d_ = m.o_cfonm_d_
    o_tcfonm_d_ = m.o_tcfonm_d_
    o_cvonm = m.o_cvonm
    o_cvonm_d_ = m.o_cvonm_d_
    o_tcvonm_d_ = m.o_tcvonm_d_
    o_cp = m.o_cp
    o_cp_d_ = m.o_cp_d_
    o_tcp_d_ = m.o_tcp_d_
    o_tpay_d_ = m.o_tpay_d_
    o_rcp_d_ = m.o_rcp_d_
    o_rcp = m.o_rcp
    o_u_d_ = m.o_u_d_
    o_u = m.o_u
    o_tu_d_ = m.o_tu_d_
    o_ehf_d_ = m.o_ehf_d_
    o_ehf = m.o_ehf
    o_tehf_d_ = m.o_tehf_d_
    o_ep0_d_ = m.o_ep0_d_
    o_ep0 = m.o_ep0
    o_tep0_d_ = m.o_tep0_d_
    o_ep1ge_d_ = m.o_ep1ge_d_
    o_ep1ge = m.o_ep1ge
    o_tep1ge_d_ = m.o_tep1ge_d_
    o_ep1gce_d_ = m.o_ep1gce_d_
    o_ep1gce = m.o_ep1gce
    o_tep1gce_d_ = m.o_tep1gce_d_
    o_ep1gcs_d_ = m.o_ep1gcs_d_
    o_ep1gcs = m.o_ep1gcs
    o_tep1gcs_d_ = m.o_tep1gcs_d_
    o_loan_last = m.o_loan_last
    o_u_cost = m.o_u_cost
    o_ehf_cost = m.o_ehf_cost
    o_ep1gcs_cost = m.o_ep1gcs_cost
    n_ep1gcs_cost = m.n_ep1gcs_cost
    r_x_in_d_ = m.r_x_in_d_
    r_x_out_d_ = m.r_x_out_d_
    r_x_o_s_d_ = m.r_x_o_s_d_
    r_x_in = m.r_x_in
    r_x_out = m.r_x_out
    o_x_in = m.o_x_in
    o_x_in_d_ = m.o_x_in_d_
    o_tx_in_d_ = m.o_tx_in_d_
    o_x_out = m.o_x_out
    o_x_out_d_ = m.o_x_out_d_
    o_tx_out_d_ = m.o_tx_out_d_
    n_x_in_d_ = m.n_x_in_d_
    n_x_out_d_ = m.n_x_out_d_
    n_x_o_s_d_ = m.n_x_o_s_d_
    n_x_in = m.n_x_in
    n_x_out = m.n_x_out
    r_ups_e_mt_in_d_ = m.r_ups_e_mt_in_d_
    r_ups_e_mt_in = m.r_ups_e_mt_in
    o_ups_e_mt_in_d_ = m.o_ups_e_mt_in_d_
    o_ups_e_mt_in = m.o_ups_e_mt_in
    o_tups_e_mt_in_d_ = m.o_tups_e_mt_in_d_
    n_ups_e_mt_in_d_ = m.n_ups_e_mt_in_d_
    n_ups_e_mt_in = m.n_ups_e_mt_in
    o_x_in_cost = m.o_x_in_cost
    n_x_in_cost = m.n_x_in_cost
    m.minimize((sum(sum(sum(p.discount[i, j] * o_pay[i, j, l] for l in L) for j in P2) for i in P) + sum(sum(sum(p.discount[i, j] * o_cfonm[i, j, l] for l in L) for j in P2) for i in P) + sum(sum(sum(p.discount[i, j] * o_cvonm[i, j, l] for l in L) for j in P2) for i in P) + sum(sum(sum(p.discount[i, j]*sum(n_pay[i0, j0, l] for i0 in P for j0 in P2 if _in((i,j),pay_window[i0, j0])) for l in L) for j in P2) for i in P) + sum(sum(sum(p.discount[i, j] * n_cfonm[i, j, l] for l in L) for j in P2) for i in P) + sum(sum(sum(p.discount[i, j] * n_cvonm[i, j, l] for l in L) for j in P2) for i in P) + sum(sum(sum(p.discount[i, j]*t_ret_cost[i, j, l] for l in L) for j in P2) for i in P) + tCoeff*sum(p.discount[n_periods,n_subperiods]* o_loan_last[n_periods, l, sT] for l in L) + tCoeff*p.discount[n_periods,n_subperiods] * sum( sum(n_loan[n_periods, n_subperiods, i0, j0, l] for j0 in P2 for i0 in P if _in((n_periods, n_subperiods),loan_window[i0, j0])  and  _ne((n_periods,n_subperiods),(i0,j0))) for l in L) + sum(p.discount[i, j] * o_u_cost[i, j, l] for l in L for j in P2 for i in P) + sum(p.discount[i, j] * n_u_cost[i, j, l] for l in L for j in P2 for i in P) + sum(p.discount[i, j] * o_ehf_cost[i, j, l] for l in L for j in P2 for i in P) + sum(p.discount[i, j] * n_ehf_cost[i, j, l] for l in L for j in P2 for i in P) + sum(p.discount[i, j] * o_ep1gcs_cost[i, j, l] for l in L for j in P2 for i in P) + sum(p.discount[i, j] * n_ep1gcs_cost[i, j, l] for l in L for j in P2 for i in P) + sum(p.discount[i, j] * o_x_in_cost[i, j, l] for l in L for j in P2 for i in P) + sum(p.discount[i, j] * n_x_in_cost[i, j, l] for l in L for j in P2 for i in P) ))
    return m

def attachPeriodBlock(m, p, s):
    P = s.P
    P2 = s.P2
    L = s.L
    Kr = s.Kr
    Kn = s.Kn
    Fu_r = s.Fu_r
    Fu_n = s.Fu_n
    Nf = s.Nf
    Nd = s.Nd
    Mt = s.Mt
    Ln = s.Ln
    n_periods = p.n_periods
    n_subperiods = p.n_subperiods
    yr_subperiod = p.yr_subperiod
    loan_subperiods = 30 // yr_subperiod
    loan_periods = loan_subperiods // n_subperiods
    remain_loan_years = 30 - loan_subperiods * yr_subperiod
    remain_subp_loan = loan_subperiods - loan_periods * n_subperiods
    remain_sub_loan = remain_subp_loan + int(remain_loan_years > 0)
    sT = 1
    sF = 2
    fP = min(P)
    fP2 = min(P2)
    fKr = min(Kr)
    fKn = min(Kn)
    r_filter = p.r_filter
    n_filter = p.n_filter
    nd_en_fltr = p.nd_en_fltr
    nd_em_fltr = p.nd_em_fltr
    node_mat = p.node_mat
    skip_mb = p.skip_mb
    links_list = p.links_list
    input_mat = p.input_mat
    output_mat = p.output_mat
    ckey = p.ckey
    key_node = p.key_node
    loan_window, pay_window = gen_lp_window(p, s)
    y_o = m.y_o
    y_r = m.y_r
    y_e = m.y_e
    e_c_d_ = m.e_c_d_
    e_c = m.e_c
    e_l = m.e_l
    x = m.x
    x_d_ = m.x_d_
    r_cp = m.r_cp
    r_cp_d_ = m.r_cp_d_
    cpb = m.cpb
    r_cpb_d_ = m.r_cpb_d_
    r_eh = m.r_eh
    r_eh_d_ = m.r_eh_d_
    r_ehf = m.r_ehf
    r_ehf_d_ = m.r_ehf_d_
    r_u = m.r_u
    r_u_d_ = m.r_u_d_
    r_u_onsite = m.r_u_onsite
    r_u_onsite_d_ = m.r_u_onsite_d_
    r_u_ehf = m.r_u_ehf
    r_u_ehf_d_ = m.r_u_ehf_d_
    r_cpe = m.r_cpe
    r_cpe_d_ = m.r_cpe_d_
    r_fu_e = m.r_fu_e
    r_fu_e_d_ = m.r_fu_e_d_
    r_u_fu_e = m.r_u_fu_e
    r_u_fu_e_d_ = m.r_u_fu_e_d_
    r_ep0 = m.r_ep0
    r_ep0_d_ = m.r_ep0_d_
    r_ep1ge = m.r_ep1ge
    r_ep1ge_d_ = m.r_ep1ge_d_
    r_ep1gce = m.r_ep1gce
    r_ep1gce_d_ = m.r_ep1gce_d_
    r_ep1gcs = m.r_ep1gcs
    r_ep1gcs_d_ = m.r_ep1gcs_d_
    r_fstck = m.r_fstck
    r_fstck_d_ = m.r_fstck_d_
    r_cfonm_d_ = m.r_cfonm_d_
    r_cfonm = m.r_cfonm
    r_cvonm_d_ = m.r_cvonm_d_
    r_cvonm = m.r_cvonm
    r_loan = m.r_loan
    r_l0 = m.r_l0
    r_l0_d_ = m.r_l0_d_
    r_l0_pd_ = m.r_l0_pd_
    r_l0add = m.r_l0add
    r_l0add_d_ = m.r_l0add_d_
    r_e_c_d_ = m.r_e_c_d_
    r_le = m.r_le
    r_le_d_ = m.r_le_d_
    r_le_pd_ = m.r_le_pd_
    r_leadd = m.r_leadd
    r_leadd_d_ = m.r_leadd_d_
    r_le_ped_ = m.r_le_ped_
    r_leadde = m.r_leadde
    r_leadde_d_ = m.r_leadde_d_
    r_pay = m.r_pay
    r_paye = m.r_paye
    e_ladd = m.e_ladd
    e_ladd_d_ = m.e_ladd_d_
    e_l_d_ = m.e_l_d_
    e_l_pd_ = m.e_l_pd_
    e_ann = m.e_ann
    e_ann_d_ = m.e_ann_d_
    e_ann_0 = m.e_ann_0
    e_ann_1 = m.e_ann_1
    e_pay = m.e_pay
    e_loan = m.e_loan
    t_loan_d_ = m.t_loan_d_
    t_ret_cost = m.t_ret_cost
    t_ret_cost_d_ = m.t_ret_cost_d_
    y_n = m.y_n
    n_c0 = m.n_c0
    n_c0_d_ = m.n_c0_d_
    n_cp = m.n_cp
    n_cp_d_ = m.n_cp_d_
    n_l_d_ = m.n_l_d_
    n_l = m.n_l
    n_ann_d_ = m.n_ann_d_
    n_ann = m.n_ann
    n_ann_0 = m.n_ann_0
    n_ann_1 = m.n_ann_1
    n_ladd_d_ = m.n_ladd_d_
    n_ladd = m.n_ladd
    n_l_pd_ = m.n_l_pd_
    n_loan = m.n_loan
    n_pay = m.n_pay
    n_eh = m.n_eh
    n_eh_d_ = m.n_eh_d_
    n_ehf = m.n_ehf
    n_ehf_d_ = m.n_ehf_d_
    n_u = m.n_u
    n_u_d_ = m.n_u_d_
    n_u_onsite = m.n_u_onsite
    n_u_onsite_d_ = m.n_u_onsite_d_
    n_u_ehf = m.n_u_ehf
    n_u_ehf_d_ = m.n_u_ehf_d_
    n_cpe = m.n_cpe
    n_cpe_d_ = m.n_cpe_d_
    n_fu_e = m.n_fu_e
    n_fu_e_d_ = m.n_fu_e_d_
    n_u_fu_e = m.n_u_fu_e
    n_u_fu_e_d_ = m.n_u_fu_e_d_
    n_ep0 = m.n_ep0
    n_ep0_d_ = m.n_ep0_d_
    n_ep1ge = m.n_ep1ge
    n_ep1ge_d_ = m.n_ep1ge_d_
    n_ep1gce = m.n_ep1gce
    n_ep1gce_d_ = m.n_ep1gce_d_
    n_ep1gcs = m.n_ep1gcs
    n_ep1gcs_d_ = m.n_ep1gcs_d_
    n_fstck = m.n_fstck
    n_fstck_d_ = m.n_fstck_d_
    n_cfonm_d_ = m.n_cfonm_d_
    n_cfonm = m.n_cfonm
    n_cvonm_d_ = m.n_cvonm_d_
    n_cvonm = m.n_cvonm
    n_u_cost = m.n_u_cost
    n_ehf_cost = m.n_ehf_cost
    o_pay = m.o_pay
    o_pay_d_ = m.o_pay_d_
    o_cfonm = m.o_cfonm
    o_cfonm_d_ = m.o_cfonm_d_
    o_tcfonm_d_ = m.o_tcfonm_d_
    o_cvonm = m.o_cvonm
    o_cvonm_d_ = m.o_cvonm_d_
    o_tcvonm_d_ = m.o_tcvonm_d_
    o_cp = m.o_cp
    o_cp_d_ = m.o_cp_d_
    o_tcp_d_ = m.o_tcp_d_
    o_tpay_d_ = m.o_tpay_d_
    o_rcp_d_ = m.o_rcp_d_
    o_rcp = m.o_rcp
    o_u_d_ = m.o_u_d_
    o_u = m.o_u
    o_tu_d_ = m.o_tu_d_
    o_ehf_d_ = m.o_ehf_d_
    o_ehf = m.o_ehf
    o_tehf_d_ = m.o_tehf_d_
    o_ep0_d_ = m.o_ep0_d_
    o_ep0 = m.o_ep0
    o_tep0_d_ = m.o_tep0_d_
    o_ep1ge_d_ = m.o_ep1ge_d_
    o_ep1ge = m.o_ep1ge
    o_tep1ge_d_ = m.o_tep1ge_d_
    o_ep1gce_d_ = m.o_ep1gce_d_
    o_ep1gce = m.o_ep1gce
    o_tep1gce_d_ = m.o_tep1gce_d_
    o_ep1gcs_d_ = m.o_ep1gcs_d_
    o_ep1gcs = m.o_ep1gcs
    o_tep1gcs_d_ = m.o_tep1gcs_d_
    o_loan_last = m.o_loan_last
    o_u_cost = m.o_u_cost
    o_ehf_cost = m.o_ehf_cost
    o_ep1gcs_cost = m.o_ep1gcs_cost
    n_ep1gcs_cost = m.n_ep1gcs_cost
    r_x_in_d_ = m.r_x_in_d_
    r_x_out_d_ = m.r_x_out_d_
    r_x_o_s_d_ = m.r_x_o_s_d_
    r_x_in = m.r_x_in
    r_x_out = m.r_x_out
    o_x_in = m.o_x_in
    o_x_in_d_ = m.o_x_in_d_
    o_tx_in_d_ = m.o_tx_in_d_
    o_x_out = m.o_x_out
    o_x_out_d_ = m.o_x_out_d_
    o_tx_out_d_ = m.o_tx_out_d_
    n_x_in_d_ = m.n_x_in_d_
    n_x_out_d_ = m.n_x_out_d_
    n_x_o_s_d_ = m.n_x_o_s_d_
    n_x_in = m.n_x_in
    n_x_out = m.n_x_out
    r_ups_e_mt_in_d_ = m.r_ups_e_mt_in_d_
    r_ups_e_mt_in = m.r_ups_e_mt_in
    o_ups_e_mt_in_d_ = m.o_ups_e_mt_in_d_
    o_ups_e_mt_in = m.o_ups_e_mt_in
    o_tups_e_mt_in_d_ = m.o_tups_e_mt_in_d_
    n_ups_e_mt_in_d_ = m.n_ups_e_mt_in_d_
    n_ups_e_mt_in = m.n_ups_e_mt_in
    o_x_in_cost = m.o_x_in_cost
    n_x_in_cost = m.n_x_in_cost
    # bm.jl:2935 — constraint o_logic_1_p_link_i_
    for i in P:
        for l in L:
            if i<n_periods:
                m.add((y_o[i, n_subperiods, l] - y_o[i+1, fP2, l]) - (0.0), '>=', 0, name='o_logic_1_p_link_i_' + ":" + repr((i, l)))
    # bm.jl:2940 — constraint r_logic_budget_p_link_i_
    for i in P:
        for l in L:
            for k in Kr:
                if k>fKr  and  i<n_periods  and  r_filter[l, k]:
                    m.add((-y_r[i, n_subperiods, l, k] + y_r[i+1, fP2, l, k]) - (0.0), '>=', 0, name='r_logic_budget_p_link_i_' + ":" + repr((i, l, k)))
    # bm.jl:2947 — constraint r_logic_onoff_1_p_link_i_
    for i in P:
        for l in L:
            for k in Kr:
                if i<n_periods  and  r_filter[l, k]:
                    m.add((y_o[i, n_subperiods, l] - y_r[i, n_subperiods, l, k] + y_r[i+1, fP2, l, k]) - (0.), '>=', 0, name='r_logic_onoff_1_p_link_i_' + ":" + repr((i, l, k)))
    # bm.jl:2953 — constraint r_logic_onoff_2_p_link_i_
    for i in P:
        for l in L:
            for k in Kr:
                if i<n_periods  and  r_filter[l, k]:
                    m.add((y_o[i, n_subperiods, l] + y_r[i, n_subperiods, l, k] - y_r[i+1, fP2, l, k]) - (0), '>=', 0, name='r_logic_onoff_2_p_link_i_' + ":" + repr((i, l, k)))
    rl0ub1d = _maximum(p.r_l0_ub, dims=2)
    # bm.jl:2963 — constraint r_l0add_m_0_p_link_i_
    for i in P:
        for l in L:
            if i>fP:
                m.add((rl0ub1d[l] * (y_r[i-1,n_subperiods,l,fKr] - y_r[i,fP2,l,fKr]) - r_l0add_d_[i, fP2, l]) - (0.), '>=', 0, name='r_l0add_m_0_p_link_i_' + ":" + repr((i, l)))
    # bm.jl:2970 — constraint r_l_m_1_y_p_link_i_
    for i in P:
        for l in L:
            if i>fP:
                m.add((rl0ub1d[l] * (-y_r[i-1, n_subperiods, l, fKr] + y_r[i,fP2,l,fKr]) - r_l0_pd_[i, fP2, l, sF]) - (-rl0ub1d[l]), '>=', 0, name='r_l_m_1_y_p_link_i_' + ":" + repr((i, l)))
    rleub1d = _maximum(p.r_le_ub, dims=2)
    # bm.jl:2976 — constraint r_leadd_m_1_p_link_i_
    for i in P:
        for l in L:
            if i>fP:
                m.add((rleub1d[l] * (y_r[i-1,n_subperiods,l,fKr] - y_r[i,fP2,l,fKr]) - r_leadd_d_[i, fP2, l]) - (0.), '>=', 0, name='r_leadd_m_1_p_link_i_' + ":" + repr((i, l)))
    # bm.jl:2981 — constraint r_le_m_1_p_link_i_
    for i in P:
        for l in L:
            if i>fP:
                m.add((rleub1d[l]*(-y_r[i-1, n_subperiods, l, fKr] + y_r[i,fP2,l,fKr]) - r_le_pd_[i, fP2, l, sF]) - (-rleub1d[l]), '>=', 0, name='r_le_m_1_p_link_i_' + ":" + repr((i, l)))
    # bm.jl:2986 — constraint r_le_m_2_p_link_i_
    for i in P:
        for l in L:
            if i>fP:
                m.add((rleub1d[l] * (y_e[i-1, n_subperiods, l] - y_e[i, fP2, l]) - r_le_ped_[i, fP2, l, sF]) - (-rleub1d[l]), '>=', 0, name='r_le_m_2_p_link_i_' + ":" + repr((i, l)))
    # bm.jl:2990 — constraint r_leadd_m_2_p_link_i_
    for i in P:
        for l in L:
            if i>fP:
                m.add((rleub1d[l] * (-y_e[i-1, n_subperiods, l] + y_e[i, fP2, l]) - r_leadde_d_[i, fP2, l]) - (0.), '>=', 0, name='r_leadd_m_2_p_link_i_' + ":" + repr((i, l)))
    # bm.jl:2996 — constraint e_logic_1_p_link_i_
    for i in P:
        for l in L:
            if i<n_periods:
                m.add((-y_e[i, n_subperiods, l] + y_e[i+1, fP2, l]) - (0), '>=', 0, name='e_logic_1_p_link_i_' + ":" + repr((i, l)))
    # bm.jl:3001 — constraint e_ladd_m_0_p_link_i_
    for i in P:
        for l in L:
            if i>fP:
                m.add((p.e_ladd_ub[l] * (-y_e[i-1, n_subperiods, l] + y_e[i, fP2, l]) - e_ladd_d_[i, fP2, l]) - (0.), '>=', 0, name='e_ladd_m_0_p_link_i_' + ":" + repr((i, l)))
    # bm.jl:3012 — constraint e_l_m_1_p_link_i_
    for i in P:
        for l in L:
            if i>fP:
                m.add((p.e_l_ub[l] * (y_e[i-1, n_subperiods, l] - y_e[i, fP2, l]) - e_l_pd_[i, fP2, l, sF]) - (-p.e_l_ub[l]), '>=', 0, name='e_l_m_1_p_link_i_' + ":" + repr((i, l)))
    # bm.jl:3030 — constraint t_ret_c_bm_0_p_link_i_
    for i in P:
        for l in L:
            if i<n_periods:
                m.add((p.t_ret_c_bM[l] * (y_o[i, n_subperiods, l] - y_o[i+1, fP2, l]) - t_ret_cost_d_[i, n_subperiods, l]) - (0.), '>=', 0, name='t_ret_c_bm_0_p_link_i_' + ":" + repr((i, l)))
    # bm.jl:3038 — constraint r_loan_d_bm_0_p_link_i_
    for i in P:
        for l in L:
            if i<n_periods:
                m.add((p.t_loan_bM[l] * (y_o[i, n_subperiods, l] - y_o[i+1, fP2, l]) - t_loan_d_[i, n_subperiods, l, sT]) - (0.), '>=', 0, name='r_loan_d_bm_0_p_link_i_' + ":" + repr((i, l)))
    # bm.jl:3045 — constraint r_loan_d_bm_1_p_link_i_
    for i in P:
        for l in L:
            if i<n_periods:
                m.add((p.t_loan_bM[l] * (-y_o[i, n_subperiods, l] + y_o[i+1, fP2, l]) - t_loan_d_[i, n_subperiods, l, sF]) - (-p.t_loan_bM[l]), '>=', 0, name='r_loan_d_bm_1_p_link_i_' + ":" + repr((i, l)))
    # bm.jl:3053 — constraint n_ladd_m_0_p_link_i_
    for i in P:
        for l in L:
            if i>fP:
                m.add((p.n_ladd_bM[l]*(y_n[i-1, n_subperiods, l, fKn] - y_n[i, fP2, l, fKn]) -n_ladd_d_[i, fP2, l]) - (0.), '>=', 0, name='n_ladd_m_0_p_link_i_' + ":" + repr((i, l)))
    # bm.jl:3059 — constraint n_ladd_m_1_p_link_i_
    for i in P:
        for l in L:
            if i>fP:
                m.add((p.n_l_bM[l] * (-y_n[i-1, n_subperiods, l, fKn] + y_n[i, fP2, l, fKn]) - n_l_pd_[i, fP2, l, sF]) - (-p.n_l_bM[l]), '>=', 0, name='n_ladd_m_1_p_link_i_' + ":" + repr((i, l)))
    # bm.jl:3066 — constraint n_logic_0_p_link_i_
    for i in P:
        for l in L:
            if i<n_periods:
                m.add((y_n[i, n_subperiods, l, fKn] - y_n[i+1, fP2, l, fKn]) - (0.), '>=', 0, name='n_logic_0_p_link_i_' + ":" + repr((i, l)))
    # bm.jl:3071 — constraint n_logic_1_p_link_i_
    for i in P:
        for l in L:
            for k in Kn:
                if i<n_periods  and  k>fKn  and  n_filter[l, k]:
                    m.add((-y_n[i, n_subperiods, l, k] + y_n[i+1, fP2, l, k]) - (0.), '>=', 0, name='n_logic_1_p_link_i_' + ":" + repr((i, l, k)))
    # bm.jl:3078 — constraint x_p_link_e_
    for i in P:
        for l in L:
            if i<n_periods:
                m.add((x[i, l] - x[i+1, l]) - (0), '==', 0, name='x_p_link_e_' + ":" + repr((i, l)))
    # bm.jl:3082 — constraint n_c0_p_link_e_
    for i in P:
        for l in L:
            if i<n_periods:
                m.add((n_c0[i, l] - n_c0[i+1, l]) - (0), '==', 0, name='n_c0_p_link_e_' + ":" + repr((i, l)))
    return m

def attachLocationBlock(m, p, s):
    P = s.P
    P2 = s.P2
    L = s.L
    Kr = s.Kr
    Kn = s.Kn
    Fu_r = s.Fu_r
    Fu_n = s.Fu_n
    Nf = s.Nf
    Nd = s.Nd
    Mt = s.Mt
    Ln = s.Ln
    n_periods = p.n_periods
    n_subperiods = p.n_subperiods
    yr_subperiod = p.yr_subperiod
    loan_subperiods = 30 // yr_subperiod
    loan_periods = loan_subperiods // n_subperiods
    remain_loan_years = 30 - loan_subperiods * yr_subperiod
    remain_subp_loan = loan_subperiods - loan_periods * n_subperiods
    remain_sub_loan = remain_subp_loan + int(remain_loan_years > 0)
    sT = 1
    sF = 2
    fP = min(P)
    fP2 = min(P2)
    fKr = min(Kr)
    fKn = min(Kn)
    r_filter = p.r_filter
    n_filter = p.n_filter
    nd_en_fltr = p.nd_en_fltr
    nd_em_fltr = p.nd_em_fltr
    node_mat = p.node_mat
    skip_mb = p.skip_mb
    links_list = p.links_list
    input_mat = p.input_mat
    output_mat = p.output_mat
    ckey = p.ckey
    key_node = p.key_node
    loan_window, pay_window = gen_lp_window(p, s)
    y_o = m.y_o
    y_r = m.y_r
    y_e = m.y_e
    e_c_d_ = m.e_c_d_
    e_c = m.e_c
    e_l = m.e_l
    x = m.x
    x_d_ = m.x_d_
    r_cp = m.r_cp
    r_cp_d_ = m.r_cp_d_
    cpb = m.cpb
    r_cpb_d_ = m.r_cpb_d_
    r_eh = m.r_eh
    r_eh_d_ = m.r_eh_d_
    r_ehf = m.r_ehf
    r_ehf_d_ = m.r_ehf_d_
    r_u = m.r_u
    r_u_d_ = m.r_u_d_
    r_u_onsite = m.r_u_onsite
    r_u_onsite_d_ = m.r_u_onsite_d_
    r_u_ehf = m.r_u_ehf
    r_u_ehf_d_ = m.r_u_ehf_d_
    r_cpe = m.r_cpe
    r_cpe_d_ = m.r_cpe_d_
    r_fu_e = m.r_fu_e
    r_fu_e_d_ = m.r_fu_e_d_
    r_u_fu_e = m.r_u_fu_e
    r_u_fu_e_d_ = m.r_u_fu_e_d_
    r_ep0 = m.r_ep0
    r_ep0_d_ = m.r_ep0_d_
    r_ep1ge = m.r_ep1ge
    r_ep1ge_d_ = m.r_ep1ge_d_
    r_ep1gce = m.r_ep1gce
    r_ep1gce_d_ = m.r_ep1gce_d_
    r_ep1gcs = m.r_ep1gcs
    r_ep1gcs_d_ = m.r_ep1gcs_d_
    r_fstck = m.r_fstck
    r_fstck_d_ = m.r_fstck_d_
    r_cfonm_d_ = m.r_cfonm_d_
    r_cfonm = m.r_cfonm
    r_cvonm_d_ = m.r_cvonm_d_
    r_cvonm = m.r_cvonm
    r_loan = m.r_loan
    r_l0 = m.r_l0
    r_l0_d_ = m.r_l0_d_
    r_l0_pd_ = m.r_l0_pd_
    r_l0add = m.r_l0add
    r_l0add_d_ = m.r_l0add_d_
    r_e_c_d_ = m.r_e_c_d_
    r_le = m.r_le
    r_le_d_ = m.r_le_d_
    r_le_pd_ = m.r_le_pd_
    r_leadd = m.r_leadd
    r_leadd_d_ = m.r_leadd_d_
    r_le_ped_ = m.r_le_ped_
    r_leadde = m.r_leadde
    r_leadde_d_ = m.r_leadde_d_
    r_pay = m.r_pay
    r_paye = m.r_paye
    e_ladd = m.e_ladd
    e_ladd_d_ = m.e_ladd_d_
    e_l_d_ = m.e_l_d_
    e_l_pd_ = m.e_l_pd_
    e_ann = m.e_ann
    e_ann_d_ = m.e_ann_d_
    e_ann_0 = m.e_ann_0
    e_ann_1 = m.e_ann_1
    e_pay = m.e_pay
    e_loan = m.e_loan
    t_loan_d_ = m.t_loan_d_
    t_ret_cost = m.t_ret_cost
    t_ret_cost_d_ = m.t_ret_cost_d_
    y_n = m.y_n
    n_c0 = m.n_c0
    n_c0_d_ = m.n_c0_d_
    n_cp = m.n_cp
    n_cp_d_ = m.n_cp_d_
    n_l_d_ = m.n_l_d_
    n_l = m.n_l
    n_ann_d_ = m.n_ann_d_
    n_ann = m.n_ann
    n_ann_0 = m.n_ann_0
    n_ann_1 = m.n_ann_1
    n_ladd_d_ = m.n_ladd_d_
    n_ladd = m.n_ladd
    n_l_pd_ = m.n_l_pd_
    n_loan = m.n_loan
    n_pay = m.n_pay
    n_eh = m.n_eh
    n_eh_d_ = m.n_eh_d_
    n_ehf = m.n_ehf
    n_ehf_d_ = m.n_ehf_d_
    n_u = m.n_u
    n_u_d_ = m.n_u_d_
    n_u_onsite = m.n_u_onsite
    n_u_onsite_d_ = m.n_u_onsite_d_
    n_u_ehf = m.n_u_ehf
    n_u_ehf_d_ = m.n_u_ehf_d_
    n_cpe = m.n_cpe
    n_cpe_d_ = m.n_cpe_d_
    n_fu_e = m.n_fu_e
    n_fu_e_d_ = m.n_fu_e_d_
    n_u_fu_e = m.n_u_fu_e
    n_u_fu_e_d_ = m.n_u_fu_e_d_
    n_ep0 = m.n_ep0
    n_ep0_d_ = m.n_ep0_d_
    n_ep1ge = m.n_ep1ge
    n_ep1ge_d_ = m.n_ep1ge_d_
    n_ep1gce = m.n_ep1gce
    n_ep1gce_d_ = m.n_ep1gce_d_
    n_ep1gcs = m.n_ep1gcs
    n_ep1gcs_d_ = m.n_ep1gcs_d_
    n_fstck = m.n_fstck
    n_fstck_d_ = m.n_fstck_d_
    n_cfonm_d_ = m.n_cfonm_d_
    n_cfonm = m.n_cfonm
    n_cvonm_d_ = m.n_cvonm_d_
    n_cvonm = m.n_cvonm
    n_u_cost = m.n_u_cost
    n_ehf_cost = m.n_ehf_cost
    o_pay = m.o_pay
    o_pay_d_ = m.o_pay_d_
    o_cfonm = m.o_cfonm
    o_cfonm_d_ = m.o_cfonm_d_
    o_tcfonm_d_ = m.o_tcfonm_d_
    o_cvonm = m.o_cvonm
    o_cvonm_d_ = m.o_cvonm_d_
    o_tcvonm_d_ = m.o_tcvonm_d_
    o_cp = m.o_cp
    o_cp_d_ = m.o_cp_d_
    o_tcp_d_ = m.o_tcp_d_
    o_tpay_d_ = m.o_tpay_d_
    o_rcp_d_ = m.o_rcp_d_
    o_rcp = m.o_rcp
    o_u_d_ = m.o_u_d_
    o_u = m.o_u
    o_tu_d_ = m.o_tu_d_
    o_ehf_d_ = m.o_ehf_d_
    o_ehf = m.o_ehf
    o_tehf_d_ = m.o_tehf_d_
    o_ep0_d_ = m.o_ep0_d_
    o_ep0 = m.o_ep0
    o_tep0_d_ = m.o_tep0_d_
    o_ep1ge_d_ = m.o_ep1ge_d_
    o_ep1ge = m.o_ep1ge
    o_tep1ge_d_ = m.o_tep1ge_d_
    o_ep1gce_d_ = m.o_ep1gce_d_
    o_ep1gce = m.o_ep1gce
    o_tep1gce_d_ = m.o_tep1gce_d_
    o_ep1gcs_d_ = m.o_ep1gcs_d_
    o_ep1gcs = m.o_ep1gcs
    o_tep1gcs_d_ = m.o_tep1gcs_d_
    o_loan_last = m.o_loan_last
    o_u_cost = m.o_u_cost
    o_ehf_cost = m.o_ehf_cost
    o_ep1gcs_cost = m.o_ep1gcs_cost
    n_ep1gcs_cost = m.n_ep1gcs_cost
    r_x_in_d_ = m.r_x_in_d_
    r_x_out_d_ = m.r_x_out_d_
    r_x_o_s_d_ = m.r_x_o_s_d_
    r_x_in = m.r_x_in
    r_x_out = m.r_x_out
    o_x_in = m.o_x_in
    o_x_in_d_ = m.o_x_in_d_
    o_tx_in_d_ = m.o_tx_in_d_
    o_x_out = m.o_x_out
    o_x_out_d_ = m.o_x_out_d_
    o_tx_out_d_ = m.o_tx_out_d_
    n_x_in_d_ = m.n_x_in_d_
    n_x_out_d_ = m.n_x_out_d_
    n_x_o_s_d_ = m.n_x_o_s_d_
    n_x_in = m.n_x_in
    n_x_out = m.n_x_out
    r_ups_e_mt_in_d_ = m.r_ups_e_mt_in_d_
    r_ups_e_mt_in = m.r_ups_e_mt_in
    o_ups_e_mt_in_d_ = m.o_ups_e_mt_in_d_
    o_ups_e_mt_in = m.o_ups_e_mt_in
    o_tups_e_mt_in_d_ = m.o_tups_e_mt_in_d_
    n_ups_e_mt_in_d_ = m.n_ups_e_mt_in_d_
    n_ups_e_mt_in = m.n_ups_e_mt_in
    o_x_in_cost = m.o_x_in_cost
    n_x_in_cost = m.n_x_in_cost
    # bm.jl:3181 — constraint ag_dem_l_link_i_
    for i in P:
        for j in P2:
            m.add((sum(o_cp[i, j, l, p.key_node] for l in L) + sum(n_cp[i, j, l, p.key_node] for l in L)) - (p.demand[i, j]), '>=', 0, name='ag_dem_l_link_i_' + ":" + repr((i, j)))
    # bm.jl:3197 — constraint ag_co2_l_link_i_
    m.add((sum( (sum(o_ep1ge[i, j, l] for l in L) + sum(o_ups_e_mt_in[i, j, l] for l in L) + sum(n_ep1ge[i, j, l] for l in L) + sum(n_ups_e_mt_in[i, j, l] for l in L) ) for j in P2 for i in P)*p.yr_subperiod) - (p.co2_budget[1, 1]), '<=', 0, name='ag_co2_l_link_i_' + ":" + repr(()))
    return m

def vintage_terms(m, p, s, init_vintage, loan0):
    P = s.P
    P2 = s.P2
    L = s.L
    Kr = s.Kr
    Kn = s.Kn
    Fu_r = s.Fu_r
    Fu_n = s.Fu_n
    Nf = s.Nf
    Nd = s.Nd
    Mt = s.Mt
    Ln = s.Ln
    n_periods = p.n_periods
    n_subperiods = p.n_subperiods
    yr_subperiod = p.yr_subperiod
    loan_subperiods = 30 // yr_subperiod
    loan_periods = loan_subperiods // n_subperiods
    remain_loan_years = 30 - loan_subperiods * yr_subperiod
    remain_subp_loan = loan_subperiods - loan_periods * n_subperiods
    remain_sub_loan = remain_subp_loan + int(remain_loan_years > 0)
    sT = 1
    sF = 2
    fP = min(P)
    fP2 = min(P2)
    fKr = min(Kr)
    fKn = min(Kn)
    r_filter = p.r_filter
    n_filter = p.n_filter
    nd_en_fltr = p.nd_en_fltr
    nd_em_fltr = p.nd_em_fltr
    node_mat = p.node_mat
    skip_mb = p.skip_mb
    links_list = p.links_list
    input_mat = p.input_mat
    output_mat = p.output_mat
    ckey = p.ckey
    key_node = p.key_node
    loan_window, pay_window = gen_lp_window(p, s)
    y_o = m.y_o
    y_r = m.y_r
    y_e = m.y_e
    e_c_d_ = m.e_c_d_
    e_c = m.e_c
    e_l = m.e_l
    x = m.x
    x_d_ = m.x_d_
    r_cp = m.r_cp
    r_cp_d_ = m.r_cp_d_
    cpb = m.cpb
    r_cpb_d_ = m.r_cpb_d_
    r_eh = m.r_eh
    r_eh_d_ = m.r_eh_d_
    r_ehf = m.r_ehf
    r_ehf_d_ = m.r_ehf_d_
    r_u = m.r_u
    r_u_d_ = m.r_u_d_
    r_u_onsite = m.r_u_onsite
    r_u_onsite_d_ = m.r_u_onsite_d_
    r_u_ehf = m.r_u_ehf
    r_u_ehf_d_ = m.r_u_ehf_d_
    r_cpe = m.r_cpe
    r_cpe_d_ = m.r_cpe_d_
    r_fu_e = m.r_fu_e
    r_fu_e_d_ = m.r_fu_e_d_
    r_u_fu_e = m.r_u_fu_e
    r_u_fu_e_d_ = m.r_u_fu_e_d_
    r_ep0 = m.r_ep0
    r_ep0_d_ = m.r_ep0_d_
    r_ep1ge = m.r_ep1ge
    r_ep1ge_d_ = m.r_ep1ge_d_
    r_ep1gce = m.r_ep1gce
    r_ep1gce_d_ = m.r_ep1gce_d_
    r_ep1gcs = m.r_ep1gcs
    r_ep1gcs_d_ = m.r_ep1gcs_d_
    r_fstck = m.r_fstck
    r_fstck_d_ = m.r_fstck_d_
    r_cfonm_d_ = m.r_cfonm_d_
    r_cfonm = m.r_cfonm
    r_cvonm_d_ = m.r_cvonm_d_
    r_cvonm = m.r_cvonm
    r_loan = m.r_loan
    r_l0 = m.r_l0
    r_l0_d_ = m.r_l0_d_
    r_l0_pd_ = m.r_l0_pd_
    r_l0add = m.r_l0add
    r_l0add_d_ = m.r_l0add_d_
    r_e_c_d_ = m.r_e_c_d_
    r_le = m.r_le
    r_le_d_ = m.r_le_d_
    r_le_pd_ = m.r_le_pd_
    r_leadd = m.r_leadd
    r_leadd_d_ = m.r_leadd_d_
    r_le_ped_ = m.r_le_ped_
    r_leadde = m.r_leadde
    r_leadde_d_ = m.r_leadde_d_
    r_pay = m.r_pay
    r_paye = m.r_paye
    e_ladd = m.e_ladd
    e_ladd_d_ = m.e_ladd_d_
    e_l_d_ = m.e_l_d_
    e_l_pd_ = m.e_l_pd_
    e_ann = m.e_ann
    e_ann_d_ = m.e_ann_d_
    e_ann_0 = m.e_ann_0
    e_ann_1 = m.e_ann_1
    e_pay = m.e_pay
    e_loan = m.e_loan
    t_loan_d_ = m.t_loan_d_
    t_ret_cost = m.t_ret_cost
    t_ret_cost_d_ = m.t_ret_cost_d_
    y_n = m.y_n
    n_c0 = m.n_c0
    n_c0_d_ = m.n_c0_d_
    n_cp = m.n_cp
    n_cp_d_ = m.n_cp_d_
    n_l_d_ = m.n_l_d_
    n_l = m.n_l
    n_ann_d_ = m.n_ann_d_
    n_ann = m.n_ann
    n_ann_0 = m.n_ann_0
    n_ann_1 = m.n_ann_1
    n_ladd_d_ = m.n_ladd_d_
    n_ladd = m.n_ladd
    n_l_pd_ = m.n_l_pd_
    n_loan = m.n_loan
    n_pay = m.n_pay
    n_eh = m.n_eh
    n_eh_d_ = m.n_eh_d_
    n_ehf = m.n_ehf
    n_ehf_d_ = m.n_ehf_d_
    n_u = m.n_u
    n_u_d_ = m.n_u_d_
    n_u_onsite = m.n_u_onsite
    n_u_onsite_d_ = m.n_u_onsite_d_
    n_u_ehf = m.n_u_ehf
    n_u_ehf_d_ = m.n_u_ehf_d_
    n_cpe = m.n_cpe
    n_cpe_d_ = m.n_cpe_d_
    n_fu_e = m.n_fu_e
    n_fu_e_d_ = m.n_fu_e_d_
    n_u_fu_e = m.n_u_fu_e
    n_u_fu_e_d_ = m.n_u_fu_e_d_
    n_ep0 = m.n_ep0
    n_ep0_d_ = m.n_ep0_d_
    n_ep1ge = m.n_ep1ge
    n_ep1ge_d_ = m.n_ep1ge_d_
    n_ep1gce = m.n_ep1gce
    n_ep1gce_d_ = m.n_ep1gce_d_
    n_ep1gcs = m.n_ep1gcs
    n_ep1gcs_d_ = m.n_ep1gcs_d_
    n_fstck = m.n_fstck
    n_fstck_d_ = m.n_fstck_d_
    n_cfonm_d_ = m.n_cfonm_d_
    n_cfonm = m.n_cfonm
    n_cvonm_d_ = m.n_cvonm_d_
    n_cvonm = m.n_cvonm
    n_u_cost = m.n_u_cost
    n_ehf_cost = m.n_ehf_cost
    o_pay = m.o_pay
    o_pay_d_ = m.o_pay_d_
    o_cfonm = m.o_cfonm
    o_cfonm_d_ = m.o_cfonm_d_
    o_tcfonm_d_ = m.o_tcfonm_d_
    o_cvonm = m.o_cvonm
    o_cvonm_d_ = m.o_cvonm_d_
    o_tcvonm_d_ = m.o_tcvonm_d_
    o_cp = m.o_cp
    o_cp_d_ = m.o_cp_d_
    o_tcp_d_ = m.o_tcp_d_
    o_tpay_d_ = m.o_tpay_d_
    o_rcp_d_ = m.o_rcp_d_
    o_rcp = m.o_rcp
    o_u_d_ = m.o_u_d_
    o_u = m.o_u
    o_tu_d_ = m.o_tu_d_
    o_ehf_d_ = m.o_ehf_d_
    o_ehf = m.o_ehf
    o_tehf_d_ = m.o_tehf_d_
    o_ep0_d_ = m.o_ep0_d_
    o_ep0 = m.o_ep0
    o_tep0_d_ = m.o_tep0_d_
    o_ep1ge_d_ = m.o_ep1ge_d_
    o_ep1ge = m.o_ep1ge
    o_tep1ge_d_ = m.o_tep1ge_d_
    o_ep1gce_d_ = m.o_ep1gce_d_
    o_ep1gce = m.o_ep1gce
    o_tep1gce_d_ = m.o_tep1gce_d_
    o_ep1gcs_d_ = m.o_ep1gcs_d_
    o_ep1gcs = m.o_ep1gcs
    o_tep1gcs_d_ = m.o_tep1gcs_d_
    o_loan_last = m.o_loan_last
    o_u_cost = m.o_u_cost
    o_ehf_cost = m.o_ehf_cost
    o_ep1gcs_cost = m.o_ep1gcs_cost
    n_ep1gcs_cost = m.n_ep1gcs_cost
    r_x_in_d_ = m.r_x_in_d_
    r_x_out_d_ = m.r_x_out_d_
    r_x_o_s_d_ = m.r_x_o_s_d_
    r_x_in = m.r_x_in
    r_x_out = m.r_x_out
    o_x_in = m.o_x_in
    o_x_in_d_ = m.o_x_in_d_
    o_tx_in_d_ = m.o_tx_in_d_
    o_x_out = m.o_x_out
    o_x_out_d_ = m.o_x_out_d_
    o_tx_out_d_ = m.o_tx_out_d_
    n_x_in_d_ = m.n_x_in_d_
    n_x_out_d_ = m.n_x_out_d_
    n_x_o_s_d_ = m.n_x_o_s_d_
    n_x_in = m.n_x_in
    n_x_out = m.n_x_out
    r_ups_e_mt_in_d_ = m.r_ups_e_mt_in_d_
    r_ups_e_mt_in = m.r_ups_e_mt_in
    o_ups_e_mt_in_d_ = m.o_ups_e_mt_in_d_
    o_ups_e_mt_in = m.o_ups_e_mt_in
    o_tups_e_mt_in_d_ = m.o_tups_e_mt_in_d_
    n_ups_e_mt_in_d_ = m.n_ups_e_mt_in_d_
    n_ups_e_mt_in = m.n_ups_e_mt_in
    o_x_in_cost = m.o_x_in_cost
    n_x_in_cost = m.n_x_in_cost
    init_vintage = dict(zip(s.L, init_vintage))
    loan0 = dict(zip(s.L, loan0))
    slots = [(i,j) for i in P for j in P2]
    init_loan_subperiods = {l: max(0, 30-init_vintage[l]) // yr_subperiod for l in L}
    init_pay_window = {l: slots[:init_loan_subperiods[l]] for l in L if init_loan_subperiods[l] > 0}
    loan0_l = {l: max(0, loan0[l] - init_vintage[l]*loan0[l]/30) for l in L}
    pay0_sp = {l: loan0_l[l]/init_loan_subperiods[l] if init_loan_subperiods[l] else 0.0 for l in L}
    S_l0 = [(i,j,l) for l in init_pay_window for i,j in init_pay_window[l]]
    S_l0_interp = [(i,j,l) for i,j,l in S_l0 if j < n_subperiods and (i,j+1) in init_pay_window[l]]
    S_l0_per = [(i,j,l) for i,j,l in S_l0 if i < n_periods and j == n_subperiods and (i+1,1) in init_pay_window[l]]
    S_loan_s_e_all = [(i,j,l) for i in P for j in P2 for l in L]
    S_loan_s_e_0 = S_l0
    S_loan_s_e = sorted(set(S_loan_s_e_all)-set(S_loan_s_e_0))
    for group in ('t_ret_c_bm_0_y_i_', 'r_loan_d_bm_0_y_i_', 'r_loan_d_bm_1_y_i_', 'r_loan_s_e_', 'o_pay_m_1_i_', 'o_tpay_s_e_', 'o_tpay_m1_i_', 'o_tpay_m0_i_'):
        m.remove_group(group)
    # bm.jl:3653 — variable l0
    l0 = m.l0 = IndexedVars()
    for _index in S_l0:
        l0[_index] = m.var('l0', _index, lb=-inf, ub=inf, binary=False)
    # bm.jl:3655 — constraint init_loanb_00
    for l in keys(init_pay_window):
        m.add((l0[(1, 1, l)]) - (loan0_l[l]), '==', 0, name='init_loanb_00' + ":" + repr((l,)))
    # bm.jl:3664 — constraint init_loanb_interperiod
    for (i,j,l) in S_l0_interp:
        m.add((l0[(i, j+1, l)]) - (l0[(i, j, l)] - pay0_sp[l]), '==', 0, name='init_loanb_interperiod' + ":" + repr((i, j, l)))
    # bm.jl:3674 — constraint init_loanb_period
    for (i,j,l) in S_l0_per:
        m.add((l0[(i+1, 1, l)]) - (l0[(i, j, l)] - pay0_sp[l]), '==', 0, name='init_loanb_period' + ":" + repr((i, j, l)))
    # bm.jl:3698 — constraint t_ret_c_bm_0_y_i_
    for i in P:
        for j in P2:
            for l in L:
                if j<n_subperiods:
                    m.add((t_ret_cost_d_[i, j, l]) - ((p.t_ret_c_bM[l] + loan0_l[l]) * (y_o[i, j, l] - y_o[i, j+1, l])), '<=', 0, name='t_ret_c_bm_0_y_i_' + ":" + repr((i, j, l)))
    # bm.jl:3706 — constraint r_loan_d_bm_0_y_i_
    for i in P:
        for j in P2:
            for l in L:
                if j<n_subperiods:
                    m.add((t_loan_d_[i, j, l, sT]) - ((p.t_loan_bM[l] + loan0_l[l])* (y_o[i, j, l] - y_o[i,j+1, l])), '<=', 0, name='r_loan_d_bm_0_y_i_' + ":" + repr((i, j, l)))
    # bm.jl:3713 — constraint r_loan_d_bm_1_y_i_
    for i in P:
        for j in P2:
            for l in L:
                if j<n_subperiods:
                    m.add((t_loan_d_[i, j, l, sF]) - ((p.t_loan_bM[l] + loan0_l[l]) * (1 + y_o[i,j+1, l] - y_o[i, j, l])), '<=', 0, name='r_loan_d_bm_1_y_i_' + ":" + repr((i, j, l)))
    # bm.jl:3725 — constraint r_loan_s_e_
    for (i,j,l) in S_loan_s_e:
        m.add((sum(r_loan[i, j, i0, j0, l] for j0 in P2 for i0 in P if _in((i, j),loan_window[i0, j0])  and  _ne((i,j),(i0,j0)) ) + sum(e_loan[i, j, i0, j0, l] for j0 in P2 for i0 in P if _in((i, j),loan_window[i0, j0])  and  _ne((i,j),(i0,j0)) ) + r_l0add[i, j, l] + r_leadd[i, j, l] + e_ladd[i, j, l]) - (t_loan_d_[i, j, l, sT] + t_loan_d_[i, j, l, sF]), '==', 0, name='r_loan_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:3737 — constraint r_loan_s_e_vint_
    for (i,j,l) in S_loan_s_e_0:
        m.add((sum(r_loan[i, j, i0, j0, l] for j0 in P2 for i0 in P if _in((i, j),loan_window[i0, j0])  and  _ne((i,j),(i0,j0)) ) + sum(e_loan[i, j, i0, j0, l] for j0 in P2 for i0 in P if _in((i, j),loan_window[i0, j0])  and  _ne((i,j),(i0,j0)) ) + r_l0add[i, j, l] + r_leadd[i, j, l] + e_ladd[i, j, l] + l0[(i,j,l)]) - (t_loan_d_[i, j, l, sT] + t_loan_d_[i, j, l, sF]), '==', 0, name='r_loan_s_e_vint_' + ":" + repr((i, j, l)))
    # bm.jl:3752 — constraint o_pay_m_1_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_pay_d_[i, j, l]) - ((p.o_pay_bM[l] + loan0_l[l]) * y_o[i, j, l]), '<=', 0, name='o_pay_m_1_i_' + ":" + repr((i, j, l)))
    # bm.jl:3758 — constraint o_tpay_s_e_
    for (i,j,l) in S_loan_s_e:
        m.add((sum(r_pay[i0, j0, l] for i0 in P for j0 in P2 if _in((i,j),pay_window[i0, j0])) + sum(r_paye[i0, j0, l] for i0 in P for j0 in P2 if _in((i,j),pay_window[i0, j0])) +sum(e_pay[i0, j0, l] for i0 in P for j0 in P2 if _in((i,j),pay_window[i0, j0]))) - (o_tpay_d_[i, j, l, sF] + o_tpay_d_[i, j, l, sT]), '==', 0, name='o_tpay_s_e_' + ":" + repr((i, j, l)))
    # bm.jl:3764 — constraint o_tpay_s_e_vint_
    for (i,j,l) in S_loan_s_e_0:
        m.add((sum(r_pay[i0, j0, l] for i0 in P for j0 in P2 if _in((i,j),pay_window[i0, j0])) + sum(r_paye[i0, j0, l] for i0 in P for j0 in P2 if _in((i,j),pay_window[i0, j0])) +sum(e_pay[i0, j0, l] for i0 in P for j0 in P2 if _in((i,j),pay_window[i0, j0])) + pay0_sp[l]) - (o_tpay_d_[i, j, l, sF] + o_tpay_d_[i, j, l, sT]), '==', 0, name='o_tpay_s_e_vint_' + ":" + repr((i, j, l)))
    # bm.jl:3775 — constraint o_tpay_m1_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_tpay_d_[i, j, l, sT]) - ((p.o_pay_bM[l] + pay0_sp[l]) * y_o[i, j, l]), '<=', 0, name='o_tpay_m1_i_' + ":" + repr((i, j, l)))
    # bm.jl:3781 — constraint o_tpay_m0_i_
    for i in P:
        for j in P2:
            for l in L:
                m.add((o_tpay_d_[i, j, l, sF]) - ((p.o_pay_bM[l] + pay0_sp[l])*(1 - y_o[i, j, l])), '<=', 0, name='o_tpay_m0_i_' + ":" + repr((i, j, l)))
    return m

def attach_supply_bounds(m, p, s):
    """Additional bounds in instance/softwareX/run_final.jl."""
    for i in s.P:
        for j in s.P2:
            captured = sum(m.o_ep1gce[i,j,l]+m.o_ep1gcs[i,j,l]+m.n_ep1gce[i,j,l]+m.n_ep1gcs[i,j,l] for l in s.L)
            m.add(captured, '<=', p.initial_co2*0.9, name=f'ccs_supply:{i},{j}')
            electricity = sum(m.o_u[i,j,l,n]+m.n_u[i,j,l,n] for l in s.L for n in s.Nd if p.nd_en_fltr[n])
            m.add(electricity*p.yr_subperiod, '<=', p.elec_budget[i,j]*p.yr_subperiod, name=f'elec_supply:{i},{j}')
    return m

def build_softwarex(p, s, init_vintage=(40,30,20,10,5), terminal_coefficient=1e-3):
    m = createBlockMod(s.P, s.L, p, s)
    attachPeriodBlock(m,p,s)
    attachLocationBlock(m,p,s)
    vintage_terms(m,p,s,init_vintage,[p.n_loanFact[l,2]*p.c0[l] for l in s.L])
    attachFullObjectiveBlock(m,p,s,terminal_coefficient)
    attach_supply_bounds(m,p,s)
    return m
