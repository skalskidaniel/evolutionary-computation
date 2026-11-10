# Lab 6 — Multiple Start Local Search (MSLS) and Iterated Local Search (ILS)
# Selective TSP with node costs (n=200, k=100).
#
# Uses the list-of-improving-moves (LM) steepest local search from Lab 5 as the
# basic local search. Implements MSLS (400 random-start LS restarts per run) and
# ILS (time budget = average MSLS running time, random-start, custom perturbation).

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

DATA = Path('data')
FIG = Path('figures')
SOL = Path('solutions')
INSTANCES = ['TSPA', 'TSPB']
RNG_SEED = 0
NUM_RUNS = 20            # runs of each method (MSLS, ILS) per instance
MSLS_RESTARTS = 400      # basic LS restarts within a single MSLS run
PERTURB_R = 5            # number of nodes removed/reinserted by the ILS perturbation

INF = np.iinfo(np.int64).max // 4


def load_instance(name: str):
    d = np.loadtxt(DATA / f'{name}.csv', delimiter=';', dtype=int)
    xy, costs = d[:, :2], d[:, 2]
    diff = xy[:, None, :] - xy[None, :, :]
    dist = np.rint(np.sqrt((diff ** 2).sum(axis=-1))).astype(np.int64)
    n = len(costs)
    k = (n + 1) // 2
    return {'name': name, 'xy': xy, 'costs': costs, 'dist': dist, 'n': n, 'k': k}


def compute_cost(tour, dist, costs) -> int:
    t = np.asarray(tour, dtype=int)
    return int(dist[t, np.roll(t, -1)].sum() + costs[t].sum())


def check_solution(tour, n: int, k: int):
    t = list(tour)
    if len(t) != k:
        return False, f'length {len(t)} != k={k}'
    if len(set(t)) != k:
        return False, 'duplicate nodes'
    if any(v < 0 or v >= n for v in t):
        return False, 'index out of range'
    return True, 'ok'


def random_tour(n: int, k: int, rng: np.random.Generator) -> list[int]:
    selected = rng.choice(n, size=k, replace=False)
    rng.shuffle(selected)
    return [int(v) for v in selected]


# ---------------------------------------------------------------------------
# List-of-improving-moves steepest local search (from Lab 5, "fair" version)
# ---------------------------------------------------------------------------

def reverse_segment(t, start, end, k):
    """Reverse t[start..end] inclusive, walking forward modulo k."""
    idxs = []
    cur = start
    while True:
        idxs.append(cur)
        if cur == end:
            break
        cur = (cur + 1) % k
    vals = [t[idx] for idx in idxs][::-1]
    for idx, v in zip(idxs, vals):
        t[idx] = v


def evaluate_all_moves_fair(t, used, D, C, k, n):
    """Evaluate the whole neighbourhood. Return improving moves (delta < 0) sorted by delta ascending."""
    moves = []
    t = np.asarray(t, dtype=int)
    pos = np.arange(k)
    prev = t[(pos - 1) % k]
    nxt = t[(pos + 1) % k]
    s = t[pos]
    U = np.flatnonzero(~used)
    m = len(U)
    if m > 0:
        inter = (D[prev[:, None], U[None, :]]
                 + D[U[None, :], nxt[:, None]]
                 - (D[prev, s] + D[s, nxt])[:, None]
                 + C[U][None, :] - C[s][:, None])
        flat = inter.ravel()
        neg = np.flatnonzero(flat < 0)
        P = np.repeat(pos, m)
        S = np.repeat(s, m)
        UU = np.tile(U, k)
        PP = np.repeat(prev, m)
        PN = np.repeat(nxt, m)
        for q in neg:
            p = int(P[q]); node_s = int(S[q]); u = int(UU[q])
            pv = int(PP[q]); nx = int(PN[q]); d = int(flat[q])
            moves.append({'kind': 'inter', 'p': p, 's': node_s, 'u': u,
                          'prev': pv, 'nxt': nx, 'delta': d,
                          'e1': (pv, node_s), 'e2': (node_s, nx)})
    ii2, jj2 = np.triu_indices(k, 2)
    valid = ~((ii2 == 0) & (jj2 == k - 1))
    ii2 = ii2[valid]; jj2 = jj2[valid]
    if len(ii2):
        a = t[ii2]; b = t[(ii2 + 1) % k]
        c = t[jj2]; dd = t[(jj2 + 1) % k]
        twoopt = D[a, c] + D[b, dd] - D[a, b] - D[c, dd]
        neg = np.flatnonzero(twoopt < 0)
        for q in neg:
            moves.append({'kind': '2opt', 'delta': int(twoopt[q]),
                          'e1': (int(a[q]), int(b[q])), 'e2': (int(c[q]), int(dd[q]))})
    moves.sort(key=lambda mv: mv['delta'])
    return moves


def move_status(mv, t, pos, k, used):
    """Return 1 (not applicable, remove from LM), 2 (different direction, keep), 3 (applicable)."""
    if mv['kind'] == 'inter':
        p = mv['p']; s = mv['s']; pv = mv['prev']; nx = mv['nxt']; u = mv['u']
        if t[p] == s and t[(p - 1) % k] == pv and t[(p + 1) % k] == nx and used[s] and not used[u]:
            return 3
        return 1
    a, b = mv['e1']; c, d = mv['e2']
    if a not in pos or b not in pos or c not in pos or d not in pos:
        return 1
    pa, pb, pc, pd = pos[a], pos[b], pos[c], pos[d]
    e1_fwd = (pb == (pa + 1) % k)
    e1_bwd = (pa == (pb + 1) % k)
    e2_fwd = (pd == (pc + 1) % k)
    e2_bwd = (pc == (pd + 1) % k)
    if e1_fwd and e2_fwd:
        if pa == pc or (pa + 1) % k == pc or pa == (pc + 1) % k:
            return 1
        return 3
    if e1_bwd and e2_bwd:
        if pb == pd or (pb + 1) % k == pd or pb == (pd + 1) % k:
            return 1
        return 3
    e1_exists = e1_fwd or e1_bwd
    e2_exists = e2_fwd or e2_bwd
    if not e1_exists or not e2_exists:
        return 1
    return 2


def apply_move(mv, t, used, k, pos):
    if mv['kind'] == 'inter':
        p = mv['p']; s = mv['s']; u = mv['u']
        used[s] = False
        used[u] = True
        t[p] = u
    else:
        a, b = mv['e1']; c, d = mv['e2']
        pa, pb = pos[a], pos[b]
        pc, pd = pos[c], pos[d]
        e1_fwd = (pb == (pa + 1) % k)
        e1_bwd = (pa == (pb + 1) % k)
        e2_fwd = (pd == (pc + 1) % k)
        e2_bwd = (pc == (pd + 1) % k)
        if e1_fwd and e2_fwd:
            start = (pa + 1) % k
            end = pc
        else:
            start = (pb + 1) % k
            end = pd
        reverse_segment(t, start, end, k)


def lm_steepest_fair(init, used, D, C, k, n):
    """Steepest local search reusing deltas via a list of improving moves (LM)."""
    t = [int(v) for v in init]
    used = used.copy()
    cost = int(compute_cost(np.array(t), D, C))
    LM = evaluate_all_moves_fair(t, used, D, C, k, n)
    n_inter = k * max(n - k, 0)
    n_2opt = k * (k - 1) // 2 - k - 1
    full_pass = n_inter + n_2opt
    n_evals = full_pass
    n_applies = 0
    while True:
        pos = {int(v): int(i) for i, v in enumerate(t)}
        applied = False
        i = 0
        while i < len(LM):
            mv = LM[i]
            st = move_status(mv, t, pos, k, used)
            if st == 1:
                LM.pop(i)
            elif st == 2:
                i += 1
            else:
                apply_move(mv, t, used, k, pos)
                cost += mv['delta']
                LM.pop(i)
                n_applies += 1
                applied = True
                break
        if applied:
            continue
        new_LM = evaluate_all_moves_fair(t, used, D, C, k, n)
        n_evals += full_pass
        if not new_LM:
            break
        LM = new_LM
    return t, int(cost), n_applies, n_evals


def used_from(tour, n):
    u = np.zeros(n, bool)
    u[np.asarray(tour, dtype=int)] = True
    return u


def steepest_ls(init, used, D, C, k, n):
    """Basic local search (LM steepest). Returns (tour, cost)."""
    tour, cost, _, _ = lm_steepest_fair(init, used, D, C, k, n)
    return tour, cost


# ---------------------------------------------------------------------------
# Perturbation for ILS: node-exchange kick
# ---------------------------------------------------------------------------

def perturb(tour, used, k, n, r, rng):
    """Remove r random tour nodes and reinsert r random unselected nodes.

    The removed positions keep their relative order for the remaining nodes;
    each removed node is replaced by a previously unselected node inserted at
    a uniformly random position. Both the subset and the order are changed,
    which kicks the search out of the current local optimum while keeping the
    kick moderate (r = 5 of k = 100 nodes) so the subsequent local search can
    still exploit the structure of the incumbent.
    """
    t = list(tour)
    used = used.copy()
    positions = sorted(rng.choice(k, size=r, replace=False), reverse=True)
    for p in positions:
        used[t[p]] = False
        del t[p]
    unselected = np.flatnonzero(~used)
    new_nodes = rng.choice(unselected, size=r, replace=False)
    for v in new_nodes:
        used[int(v)] = True
        insert_pos = int(rng.integers(0, len(t) + 1))
        t.insert(insert_pos, int(v))
    return t, used


# ---------------------------------------------------------------------------
# MSLS
# ---------------------------------------------------------------------------

def mls_run(D, C, k, n, rng):
    """Single MSLS run: MSLS_RESTARTS basic LS from random starts, keep best."""
    best_cost = np.inf
    best_tour = None
    for _ in range(MSLS_RESTARTS):
        init = random_tour(n, k, rng)
        used = np.zeros(n, bool)
        used[np.array(init)] = True
        tour, cost = steepest_ls(init, used.copy(), D, C, k, n)
        if cost < best_cost:
            best_cost = cost
            best_tour = tour.copy()
    return best_tour, best_cost


# ---------------------------------------------------------------------------
# ILS
# ---------------------------------------------------------------------------

def ils_run(D, C, k, n, time_budget, rng, r=PERTURB_R):
    """Single ILS run: random start, LS, then perturb+LS until time_budget.

    Returns (best_tour, best_cost, num_ls_runs).
    Acceptance: accept the perturbed+locally-searched solution if it is
    strictly better than the current incumbent (standard hill-climbing ILS).
    """
    init = random_tour(n, k, rng)
    used = used_from(init, n)
    tour, cost = steepest_ls(init, used.copy(), D, C, k, n)
    used = used_from(tour, n)
    best_tour = tour.copy()
    best_cost = cost
    num_ls_runs = 1
    t0 = time.perf_counter()
    while time.perf_counter() - t0 < time_budget:
        ptour, pused = perturb(tour, used, k, n, r, rng)
        qtour, qcost = steepest_ls(ptour, pused.copy(), D, C, k, n)
        num_ls_runs += 1
        if qcost < cost:
            tour, cost = qtour, qcost
            used = used_from(qtour, n)
            if cost < best_cost:
                best_cost = cost
                best_tour = tour.copy()
    return best_tour, best_cost, num_ls_runs


def mlabel(key):
    if key == 'msls':
        return 'MSLS (400 restarts/run)'
    if key == 'ils':
        return 'ILS (time budget = avg MSLS time)'
    return key


if __name__ == '__main__':
    FIG.mkdir(exist_ok=True)
    SOL.mkdir(exist_ok=True)

    problems = {name: load_instance(name) for name in INSTANCES}
    print('loaded instances:', {n: (problems[n]['n'], problems[n]['k']) for n in INSTANCES})

    results = {}
    PLOT_METHODS = ['msls', 'ils']
    best_tours = {name: {key: None for key in PLOT_METHODS} for name in INSTANCES}
    best_costs = {name: {key: None for key in PLOT_METHODS} for name in INSTANCES}

    for name in INSTANCES:
        pr = problems[name]
        n, k, D, C = pr['n'], pr['k'], pr['dist'], pr['costs']
        rng = np.random.default_rng(RNG_SEED)

        results[name] = {'n': n, 'k': k, 'methods': {}}

        # ---------------- MSLS ----------------
        t0 = time.perf_counter()
        m_costs, m_tours = [], []
        for _ in range(NUM_RUNS):
            tour, cost = mls_run(D, C, k, n, rng)
            m_costs.append(cost)
            m_tours.append(tour)
        m_time = time.perf_counter() - t0
        m_costs = np.array(m_costs)
        mj = int(np.argmin(m_costs))
        avg_msls_time = m_time / NUM_RUNS
        results[name]['methods']['msls'] = {
            'costs': m_costs, 'mean': float(m_costs.mean()),
            'min': int(m_costs.min()), 'max': int(m_costs.max()),
            'time': m_time, 'time_per_run': avg_msls_time,
            'best_tour': m_tours[mj], 'best_cost': int(m_costs.min()),
        }
        print(f"{name} MSLS: mean={m_costs.mean():.1f} ({m_costs.min()}-{m_costs.max()}) "
              f"total={m_time:.2f}s per_run={avg_msls_time:.2f}s")

        # ---------------- ILS ----------------
        time_budget = avg_msls_time
        t0 = time.perf_counter()
        i_costs, i_tours, i_nls = [], [], []
        for _ in range(NUM_RUNS):
            tour, cost, nls = ils_run(D, C, k, n, time_budget, rng)
            i_costs.append(cost)
            i_tours.append(tour)
            i_nls.append(nls)
        i_time = time.perf_counter() - t0
        i_costs = np.array(i_costs)
        ij = int(np.argmin(i_costs))
        results[name]['methods']['ils'] = {
            'costs': i_costs, 'mean': float(i_costs.mean()),
            'min': int(i_costs.min()), 'max': int(i_costs.max()),
            'time': i_time, 'time_per_run': float(i_time / NUM_RUNS),
            'budget': float(time_budget),
            'nls_runs': i_nls, 'nls_mean': float(np.mean(i_nls)),
            'nls_min': int(np.min(i_nls)), 'nls_max': int(np.max(i_nls)),
            'best_tour': i_tours[ij], 'best_cost': int(i_costs.min()),
        }
        print(f"{name} ILS: mean={i_costs.mean():.1f} ({i_costs.min()}-{i_costs.max()}) "
              f"total={i_time:.2f}s budget={time_budget:.2f}s "
              f"LSruns/ILSrun={np.mean(i_nls):.1f} ({np.min(i_nls)}-{np.max(i_nls)})")

        best_tours[name]['msls'] = results[name]['methods']['msls']['best_tour']
        best_costs[name]['msls'] = results[name]['methods']['msls']['best_cost']
        best_tours[name]['ils'] = results[name]['methods']['ils']['best_tour']
        best_costs[name]['ils'] = results[name]['methods']['ils']['best_cost']

    # ---------------- Results tables ----------------
    def mlabel2(key):
        if key == 'msls':
            return 'MSLS (400 restarts/run)'
        return 'ILS (time budget = avg MSLS time)'

    rows_obj, rows_time, rows_ils = [], [], []
    for name in INSTANCES:
        for key in ['msls', 'ils']:
            m = results[name]['methods'][key]
            rows_obj.append({'Instance': name, 'Method': mlabel2(key),
                             'mean': m['mean'], 'min': m['min'], 'max': m['max']})
            rows_time.append({'Instance': name, 'Method': mlabel2(key),
                              'total [s]': m['time'], 'per run [s]': m['time_per_run']})
            if key == 'ils':
                rows_ils.append({'Instance': name, 'Method': mlabel2(key),
                                 'budget [s]': m['budget'],
                                 'LS runs/ILS run (mean)': m['nls_mean'],
                                 'LS runs/ILS run (min)': m['nls_min'],
                                 'LS runs/ILS run (max)': m['nls_max']})

    print()
    print('Objective function: mean (min-max), 20 runs per instance')
    print(pd.DataFrame(rows_obj).round({'mean': 1}).to_string(index=False))
    print()
    print('Running time, 20 runs per instance')
    print(pd.DataFrame(rows_time).round({'total [s]': 2, 'per run [s]': 3}).to_string(index=False))
    print()
    print('ILS: number of runs of basic local search')
    print(pd.DataFrame(rows_ils).round({'budget [s]': 2, 'LS runs/ILS run (mean)': 1}).to_string(index=False))

    # ---------------- Reference: best local search from Lab 5 ----------------
    # Values taken from Lab 5, 200 runs unless noted.
    ref_rows = [
        {'Instance': 'TSPA', 'Method': 'Best local search (lm_steepest, Lab 5, ref)',
         'mean': 76679.1, 'min': 73284, 'max': 80508},
        {'Instance': 'TSPA', 'Method': 'Full steepest 2-opt (baseline, Lab 5, ref)',
         'mean': 76278.1, 'min': 73843, 'max': 79671},
        {'Instance': 'TSPB', 'Method': 'Best local search (lm_steepest, Lab 5, ref)',
         'mean': 52708.5, 'min': 50202, 'max': 55616},
        {'Instance': 'TSPB', 'Method': 'Full steepest 2-opt (baseline, Lab 5, ref)',
         'mean': 52598.8, 'min': 50586, 'max': 55376},
    ]
    print()
    print('Reference values from previous labs (200 runs, not part of the 20-run experiment):')
    print(pd.DataFrame(ref_rows).round({'mean': 1}).to_string(index=False))

    # ---------------- Visualization ----------------
    all_costs = np.concatenate([problems[n]['costs'] for n in INSTANCES])
    norm = Normalize(vmin=int(all_costs.min()), vmax=int(all_costs.max()))
    sm = plt.cm.ScalarMappable(norm=norm, cmap='viridis')

    for name in INSTANCES:
        for key in PLOT_METHODS:
            pr = problems[name]
            xy, costs = pr['xy'], pr['costs']
            tour = best_tours[name][key]
            cost = best_costs[name][key]
            fig, ax = plt.subplots(figsize=(6.0, 4.5))
            ax.set_facecolor('white')
            ax.scatter(xy[:, 0], xy[:, 1], c=costs, cmap='viridis', norm=norm, s=14, alpha=0.9, zorder=1)
            path = xy[np.append(tour, tour[0])]
            ax.plot(path[:, 0], path[:, 1], 'k-', lw=0.8, zorder=3)
            ax.scatter(path[:, 0], path[:, 1], c=costs[np.append(tour, tour[0])],
                       cmap='viridis', norm=norm, s=22, edgecolor='black', linewidth=0.3, zorder=4)
            ax.set_title(f'{name} / {mlabel2(key)}\nbest={cost}', fontsize=11)
            ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
            fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.04, label='node cost')
            fname = FIG / f'lab6_{name}_{key}.png'
            fig.savefig(fname, dpi=150, bbox_inches='tight')
            plt.close(fig)
            print(f'saved {fname}')

    # ---------------- Best solutions as node index lists ----------------
    for name in INSTANCES:
        print(f'=== {name} (n={problems[name]["n"]}, k={problems[name]["k"]}) ===')
        for key in PLOT_METHODS:
            tour = best_tours[name][key]
            cost = best_costs[name][key]
            print(f'{mlabel2(key)}: cost={cost}')
            print(tour)
            print()

    # ---------------- Save solution CSVs ----------------
    for name in INSTANCES:
        for key in PLOT_METHODS:
            tour = [int(v) for v in best_tours[name][key]]
            with open(SOL / f'lab6_{name}_{key}_ids.csv', 'w') as f:
                for idx in tour:
                    f.write(f'{idx}\n')
            print(f'saved solutions/lab6_{name}_{key}_ids.csv')

    # ---------------- Solution checking ----------------
    print()
    for name in INSTANCES:
        pr = problems[name]
        n, k = pr['n'], pr['k']
        for key in PLOT_METHODS:
            tour = best_tours[name][key]
            ok, msg = check_solution(tour, n, k)
            cost = compute_cost(tour, pr['dist'], pr['costs'])
            stored = best_costs[name][key]
            print(f'{name}/{key}: ok={ok} ({msg}) recomputed_cost={cost} stored={stored} match={cost == stored}')
