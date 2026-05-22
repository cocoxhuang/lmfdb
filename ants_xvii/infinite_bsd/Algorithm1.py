"""Algorithm1.py

This script generates a text file which includes elliptic curves with
a family of quadratic twists that satisfy the full BSD conjecture formula.
It corresponds to Algorithm 1 in the paper.

To reproduce the paper's full output file (curves over the entire LMFDB,
with the BSD(E,2) verification step active), run the following command
from this level of the directory in the terminal:

    PYTHONUNBUFFERED=1 sage -python -u Algorithm1.py --ec_file output/ec_labels_500k.txt

PYTHONUNBUFFERED=1 (and -u) ensure progress prints land in the log without
buffering; expected runtime is roughly 8 minutes. The output is written to
the path passed via --ec_file; a sibling <base>_fail<ext> file is written
only if any curve fails the BSD(E,2) check.

To limit the search to curves with conductor < 150 (recommended for testing):

    sage -python Algorithm1.py --cond_upper_bound 150

"""

import sys
import os
# Add repo root to path (go up 2 levels: infinite_bsd -> ants_xvii -> repo_root)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import argparse
from lmfdb import db
import pandas as pd
import numpy as np
import pytz
from sage.all import (EllipticCurve, primes_first_n, round, Integer, QQ, RR,
                    ZZ, primes, fundamental_discriminant, gcd, kronecker_symbol,
                    PolynomialRing)
from sage.schemes.elliptic_curves.BSD import (
    mwrank_two_descent_work,
    pari_two_descent_work,
    native_two_isogeny_descent_work,
)
import multiprocessing as mp
import time
from datetime import datetime

MWRANK_DEFAULT_TIMEOUT_SECS = 30

# Load the elliptic curve table from the LMFDB
ecq = db.ec_curvedata
ec_classdata = db.ec_classdata
ec_mwbsd = db.ec_mwbsd

# Information columns from lmfdb
ECQ_COLS = ['ainvs', 'lmfdb_label', 'conductor', 'rank', 'torsion', 'absD', 'bad_primes',
'regulator', 'sha', 'lmfdb_iso', 'torsion_structure', 'torsion_primes', 'Clabel', 'isogeny_degrees']
MWBSD_COLS = ['lmfdb_label', 'tamagawa_product'] + ['real_period', 'special_value']

def get_current_time_str():
    # Get the current time in UTC
    utc_now = datetime.now(pytz.utc)

    # Convert UTC time to US Eastern Time
    eastern = pytz.timezone('US/Eastern')
    eastern_now = utc_now.astimezone(eastern)

    # Format the time in a human-readable format
    return eastern_now.strftime("Generated at %H:%M (eastern) on %A %d %B %Y")

def is_ramified(bad_primes: list[int], minimal_disc: Integer) -> bool:
    '''returns True if the curve is ramified at any bad primes'''
    return all(
        any(minimal_disc.valuation(p) % ell != 0 for p in bad_primes if p != ell)
        for ell in bad_primes
    )

def merge_a3(df: pd.DataFrame) -> pd.DataFrame:
    '''merges a3 (=a_p at p=3) from ec_classdata into the given dataframe.'''
    payload = ec_classdata.search({'lmfdb_iso': {'$in': df['lmfdb_iso'].tolist()}},
                                  projection=['lmfdb_iso', 'aplist'])
    cd = pd.DataFrame(list(payload))
    cd['a3'] = cd['aplist'].apply(lambda x: x[1])
    return pd.merge(df, cd[['lmfdb_iso', 'a3']], on='lmfdb_iso', how='inner')

def merge_mwbsd(df: pd.DataFrame) -> pd.DataFrame:
    '''merges mwbsd data into the given dataframe'''
    mwbsd_query = {'lmfdb_label': {'$in': df['lmfdb_label'].tolist()}}
    mwbsd_payload = ec_mwbsd.search(mwbsd_query, projection=MWBSD_COLS)
    mwbsd_df = pd.DataFrame(list(mwbsd_payload))
    df = pd.merge(df, mwbsd_df, on='lmfdb_label', how='inner')
    return df

def filter_CONDITION_p_isogeny(df: pd.DataFrame) -> pd.DataFrame:
    '''CONDITION: for every p in A, E admits no rational p-isogeny, where
        A := {3 if 3 | N or |a3| = 3} \\cup {5 if 5 | N} \\cup {7 if 7 | N}.
    N is the conductor (equivalently, p | N iff p in bad_primes).'''
    data_idx = []
    for index, curve in df.iterrows():
        bad_primes = set(curve['bad_primes'])
        A = set()
        if 3 in bad_primes or abs(curve['a3']) == 3:
            A.add(3)
        if 5 in bad_primes:
            A.add(5)
        if 7 in bad_primes:
            A.add(7)
        isogeny_primes = [p for p in curve['isogeny_degrees'] if Integer(p).is_prime()]
        if all(p not in isogeny_primes for p in A):
            data_idx.append(index)
    return df.loc[data_idx]

def filter_CONDITION_ramification(df: pd.DataFrame) -> pd.DataFrame:
    '''CONDITION: ramification at any bad primes'''
    data_idx = []
    for index,curve in df.iterrows():
        ainvs = curve['ainvs']
        minimal_disc = Integer(curve['absD'])
        bad_primes = curve['bad_primes']

        CONDITION = is_ramified(bad_primes, minimal_disc)

        if CONDITION:
            data_idx.append(index)
    return df.loc[data_idx]

def filter_CONDITION_2_torsion(df: pd.DataFrame) -> pd.DataFrame:
    '''filters dataframe for CONDITION: E(Q)[2] = Z/2Z'''
    allowed_torsion_structures = [[2], [4], [6], [8], [10], [12]]
    df = df[df['torsion_structure'].isin(allowed_torsion_structures)] 
    data_idx = []
    for index,curve in df.iterrows():
        ainvs = curve['ainvs']
        E = EllipticCurve(ainvs)
        E_torsion_gens = E.torsion_subgroup().gens()
        if len(E_torsion_gens) == 1:
            data_idx.append(index)
    return df.loc[data_idx]

def filter_CONDITION_E_prime(df: pd.DataFrame) -> pd.DataFrame:
    '''CONDITION: Sha(E')[2] = 0 and E'(Q)[2] = Z/2Z, where E' = E / E(Q)[2].

    Sha(E')[2] = 0 is verified unconditionally via 2-descent on E':
    check_BSD_at_2(E', 1, sha_an_ord_2(E')) returns True iff
    sha_an_ord_2(E') = 0 AND descent on E' pins dim_{F_2} Sha(E')[2] = 0,
    i.e., Sha(E')[2] = 0 algebraically. The analytic value sha_an(E') from
    E_prime.sha().an() is used only as input to that check; it is not the
    statement being asserted.
    '''
    data_idx = []
    for index, curve in df.iterrows():
        ainvs = curve['ainvs']
        E = EllipticCurve(ainvs)

        # The unique order-2 generator of E(Q)_tors: (n/2) * P where P generates
        # the cyclic torsion of order n. Works uniformly for n in {2,4,6,8,10,12}.
        E_torsion_gens = E.torsion_subgroup().gens()
        assert len(E_torsion_gens) == 1
        P = E_torsion_gens[0]
        order = P.order()
        if order % 2 != 0:
            continue
        E_two_torsion_gen = (order // 2) * P
        assert E_two_torsion_gen.order() == 2

        # E' = E / <E(Q)[2]>
        C = E(E_two_torsion_gen)
        E_prime = E.isogeny_codomain(C)

        # First condition: E'(Q)[2] = Z/2Z, i.e. E'(Q)_tors cyclic of even order.
        n_gen = E_prime.torsion_subgroup().ngens()
        E_prime_tors_order = E_prime.torsion_order()
        if not (n_gen == 1 and E_prime_tors_order % 2 == 0):
            continue

        # Analytic Sha order for E', used as the sha_an_ord_2 input below.
        try:
            E_prime_sha_an = E_prime.sha().an().round()
        except Exception as e:
            print(f"Error calculating sha(E') for {curve['lmfdb_label']}: {e}")
            try:
                E_prime_sha_an = ecq.lucky({'ainvs': [int(x) for x in E_prime.ainvs()]},
                                           projection='sha')
            except Exception:
                continue  # cannot get sha_an for E'; reject this row
        sha_an_ord_2_E_prime = ZZ(E_prime_sha_an).valuation(2)

        # Second condition: Sha(E')[2] = 0, verified rigorously via 2-descent on E'.
        # check_BSD_at_2 returns True iff sha_an_ord_2 = 0 AND descent pins dim Sha[2] = 0.
        if check_BSD_at_2(E_prime, two_tor_rk=1, sha_an_ord_2=sha_an_ord_2_E_prime):
            data_idx.append(index)
    return df.loc[data_idx]

def filter_CONDITION_S(df: pd.DataFrame, S_bound: int = 10000) -> pd.DataFrame:
    '''CONDITION: S is nonempty where S = { q = 1 mod 4 : q does not divide N, ord_2(N_q)) = 1 }'''
    data_idx = []
    for index,curve in df.iterrows():
        ainvs = curve['ainvs']
        E = EllipticCurve(ainvs)
        cond = curve['conductor']

        candidates = list(primes(S_bound))
        for q in candidates:
            if q % 4 == 1 and cond % q != 0:
                N_q = ZZ(1 + q - E.ap(q))
                if N_q.valuation(2) == 1:
                    data_idx.append(index)
                    break
    return df.loc[data_idx]

def filter_CONDITION_S_deterministic(df: pd.DataFrame) -> pd.DataFrame:
    '''Deterministic counterpart of filter_CONDITION_S.

    Assumes the input rows have E(Q)[2] = Z/2Z (enforced upstream by
    filter_CONDITION_2_torsion), so f(x) = x^3 + A x + B (from the short
    Weierstrass model y^2 = x^3 + A x + B) has exactly one rational root x0.
    S is empty iff any of {3*x0^2 + A, -(3*x0^2 + A), 4*A^3 + 27*B^2} is a
    rational square; otherwise S is nonempty and the row is kept.
    '''
    P = PolynomialRing(QQ, 'x')
    data_idx = []
    for index, curve in df.iterrows():
        ainvs = curve['ainvs']
        E = EllipticCurve(ainvs).short_weierstrass_model()
        A = E.a4()
        B = E.a6()
        f = P([B, A, 0, 1])
        roots = f.roots()
        assert len(roots) == 1, (
            f"Expected exactly one rational root for {curve['lmfdb_label']}, got {len(roots)}"
        )
        x0 = roots[0][0]
        num1 = 3 * x0**2 + A
        num2 = 4 * A**3 + 27 * B**2
        if not (num1.is_square() or (-num1).is_square() or num2.is_square()):
            data_idx.append(index)
    return df.loc[data_idx]

def ord_2_two_torsion_order(torsion_structure, torsion):
    '''returns ord_2( the order of the 2-torsion subgroup)'''
    if len(torsion_structure) == 1:
        two_torsion_val_2 = min(ZZ(torsion).valuation(2), 1)
    else:
        two_torsion_val_2 = ZZ(torsion_structure[0]).valuation(2) + ZZ(torsion_structure[1]).valuation(2)
    return two_torsion_val_2

def two_torsion_rank_from_structure(torsion_structure) -> int:
    '''dim_{F_2} E(Q)[2] from the LMFDB torsion_structure field.

    Possible torsion_structure values for E/Q: [] or [1] (trivial), [n] (cyclic),
    or [m, n] with m | n (necessarily Z/2 x Z/2k by Mazur).
    '''
    if not torsion_structure or torsion_structure == [1]:
        return 0
    if len(torsion_structure) == 1:
        return 1 if torsion_structure[0] % 2 == 0 else 0
    return 2

def _mwrank_descent_worker(ainvs, two_tor_rk, queue):
    '''Subprocess entry point: runs mwrank's 2-descent and pushes the
    (sha2_lo, sha2_hi) bounds onto the queue. Used by
    mwrank_two_descent_with_timeout to give mwrank a hard wall-clock budget.'''
    try:
        E = EllipticCurve(list(ainvs))
        _, _, sha2_lo, sha2_hi, _ = mwrank_two_descent_work(E, two_tor_rk)
        queue.put(('ok', int(sha2_lo), int(sha2_hi)))
    except Exception as e:
        queue.put(('err', repr(e)))

def mwrank_two_descent_with_timeout(E, two_tor_rk: int, timeout_secs: float):
    '''Run mwrank's 2-descent in a forked subprocess; raise TimeoutError on overrun.

    mwrank's C++ inner loop (rank1::getquartics) doesn't yield to Python, so
    signal-based timeouts can't interrupt it. We fork a worker, join with a
    deadline, and SIGTERM (then SIGKILL) on overrun.
    '''
    ctx = mp.get_context('fork')
    queue = ctx.Queue()
    ainvs = tuple(int(a) for a in E.ainvs())
    proc = ctx.Process(target=_mwrank_descent_worker, args=(ainvs, two_tor_rk, queue))
    proc.start()
    proc.join(timeout=timeout_secs)
    if proc.is_alive():
        proc.terminate()
        proc.join(timeout=2)
        if proc.is_alive():
            proc.kill()
            proc.join(timeout=2)
        raise TimeoutError(f"mwrank exceeded {timeout_secs}s for ainvs={ainvs}")
    try:
        msg = queue.get_nowait()
    except Exception:
        raise RuntimeError("mwrank worker exited without result")
    if msg[0] == 'ok':
        return msg[1], msg[2]
    raise RuntimeError(f"mwrank worker error: {msg[1]}")

def check_BSD_at_2(E, two_tor_rk: int, sha_an_ord_2: int,
                   mwrank_timeout_secs: float = MWRANK_DEFAULT_TIMEOUT_SECS) -> bool:
    '''Return True iff the 2-part of BSD for E is unconditionally proved by 2-descent.

    BSD(E,2) is the equality ord_2(#Sha(E)) = ord_2(#Sha_an(E)). 2-descent
    computes dim_{F_2} Sha(E)[2], not ord_2(#Sha(E)) directly. By Cassels-Tate
    we only have dim Sha[2] <= ord_2(#Sha), with equality iff Sha[2^infinity]
    is elementary abelian. So 2-descent alone settles BSD(E,2) only when
    sha_an_ord_2 = 0: a descent that pins dim Sha[2] = 0 forces Sha[2] = 0,
    hence Sha[2^infinity] = 0, hence ord_2(#Sha) = 0 = sha_an_ord_2.

    For sha_an_ord_2 >= 1, matching descent bounds to the analytic value would
    only verify a necessary condition (dim Sha[2] = ord_2(#Sha_an)); it would
    NOT rule out a non-elementary Sha[2^infinity] (e.g. (Z/4)^2 with
    dim Sha[2] = 2 but ord_2(#Sha) = 4 != ord_2(#Sha_an) = 2). Settling
    BSD(E,2) in that regime requires a higher descent (4-descent, 8-descent,
    ...), which this function does not perform. To avoid overclaiming, the
    function therefore returns False whenever sha_an_ord_2 != 0, irrespective
    of what 2-descent reports. Sage's prove_BSD prints "p = 2: True by
    2-descent" in this regime as well; that claim is stronger than the descent
    actually establishes, and this function does not follow it.

    When sha_an_ord_2 = 0, cascades through pari -> mwrank -> native sage
    two-isogeny descent; returns True as soon as any backend pins both bounds
    to 0. Pari runs first because its 2-descent is fast and resolves the 2-part
    for the great majority of curves in our pipeline; mwrank is strictly more
    powerful but slow and is the only backend that can wedge for many minutes
    on pathological inputs, so it runs in a forked subprocess with a wall-clock
    timeout (mwrank_timeout_secs).

    Assumes E is already the optimal curve in its isogeny class and that the
    caller supplies two_tor_rk = dim_{F_2} E(Q)[2] and sha_an_ord_2 = ord_2 of
    the analytic Sha (e.g. pulled from the LMFDB), so this function does no
    redundant lookups.
    '''
    if sha_an_ord_2 != 0:
        return False
    descent_calls = (
        lambda: pari_two_descent_work(E)[2:4],
        lambda: mwrank_two_descent_with_timeout(E, two_tor_rk, mwrank_timeout_secs),
        lambda: native_two_isogeny_descent_work(E, two_tor_rk)[2:4],
    )
    for call in descent_calls:
        try:
            sha2_lo, sha2_hi = call()
        except Exception:
            continue
        if sha2_lo == 0 == sha2_hi:
            return True
    return False

# Main function
def bsd_infinite_twists(cond_upper_bound:int = None, cond_lower_bound:int = None,
                        EC_FILE = 'output/ec_labels.txt',
                        skip_filter_S = False,
                        use_deterministic_S = True,
                        skip_BSD_at_2_check = False):
    '''
    Generates two text files containing elliptic curves over Q
    which have an infinite family of quadratic twists satisfying the full BSD conjecture formula.
    The first file contains LMFDB labels, and the second file contains Cremona labels.

    2-part of the BSD verification references:
        [CLZ20] : L. Cai, C. Li and S. Zhai, On the 2-part of the Birch and Swinnerton-Dyer conjecture for quadratic twists of elliptic
    curves, J. Lond. Math. Soc. (2) 101 (2020), no. 2, 714–734.
        [Zha16] :  S. Zhai, Non-vanishing theorems for quadratic twists of elliptic curves, Asian J. Math. 20 (2016), no. 3, 475–502
    
    Args:
        cond_upper_bound (int): Upper bound for the conductor of elliptic curves to consider. Defaults to None.
        cond_lower_bound (int): Lower bound for the conductor of elliptic curves to consider. Defaults to None.
        EC_FILE (str): Path to save the elliptic curve labels file. Defaults to 'output/ec_labels.txt'
        skip_filter_S (bool): FOR TESTING ONLY. If True, skips the S-nonempty filter on the CLZ20 branch.
            This is unsound: skipping S lets through curves that should be filtered out. Defaults to False.
            Do not enable for production runs.
        use_deterministic_S (bool): If True (default), uses the deterministic S filter
            (filter_CONDITION_S_deterministic) — the form stated in the paper's algorithm.
            If False, uses the probabilistic filter_CONDITION_S. Has no effect if
            skip_filter_S is True.
        skip_BSD_at_2_check (bool): If True, skips the unconditional BSD(E,2) verification step
            (check_BSD_at_2). Defaults to False (the check runs). Intended for faster re-runs when
            BSD(E,2) has already been verified on the same curves in a prior run; with this flag set,
            output curves are only guaranteed to satisfy the upstream filters and are NOT verified
            for BSD(E,2). No fail file is written in this mode.
    '''
    start_time = time.time()
    run_summary = f"Params: cond_upper_bound={cond_upper_bound}, cond_lower_bound={cond_lower_bound}"
    
    # Path to save the output files
    ec_file = EC_FILE

    # Get curve data from the LMFDB
    ecq_query = {
                'semistable' : True,  # CONDITION 1a: semistable
                'num_bad_primes' : {'$gte' : 2},  # CONDITION 1d: at least two bad prime
                'optimality' : 1,  # CONDITION 1e: E is optimal
                }
    if cond_upper_bound is not None:
        if cond_lower_bound is not None:
            ecq_query['conductor'] = {'$gte' : cond_lower_bound,  '$lt': cond_upper_bound}
        else:
            ecq_query['conductor'] = {'$lt': cond_upper_bound}
    else:
        if cond_lower_bound is not None:
            ecq_query['conductor'] = {'$gte' : cond_lower_bound}
    ecq_payload = ecq.search(ecq_query, projection=ECQ_COLS)
    df = pd.DataFrame(list(ecq_payload))
    assert df['lmfdb_iso'].nunique() == len(df), "Values in the 'lmfdb_iso' column are not unique!"

    # pull a3 from ec_classdata; it gates whether 3 enters the isogeny-check set A
    df = merge_a3(df)

    # CONDITION: no rational p-isogeny for p in A (see filter_CONDITION_p_isogeny)
    df = filter_CONDITION_p_isogeny(df)

    # CONDITION: ramification at ANY bad primes
    df = filter_CONDITION_ramification(df)

    # checking the 2-part of BSD
    df = df[df['rank'] == 0]    # only rank 0 curves
    # --------------------------------------------------------------
    # using criteria from [CLZ20, Theorem 1.5]
    # --------------------------------------------------------------

    # CONDITION: E(Q)[2] = Z/2Z
    df_CLZ20 = filter_CONDITION_2_torsion(df)
    
    # CONDITION: 2-ord of special_value/(real_period*regulator) = -1
    df_CLZ20 = merge_mwbsd(df_CLZ20)
    df_CLZ20['L_alg'] = (df_CLZ20['special_value']/(df_CLZ20['real_period'] * df_CLZ20['regulator']))
    df_CLZ20['L_alg'] = df_CLZ20['L_alg'].apply(lambda x: QQ(RR(x)).valuation(2))    # now L_alg is 2-valuation of the quantity
    df_CLZ20 = df_CLZ20[df_CLZ20['L_alg'] == -1]  
    df_CLZ20 = filter_CONDITION_E_prime(df_CLZ20)

    # CONDITION: S is nonempty (up to a bound)
    if not skip_filter_S:
        if use_deterministic_S:
            df_CLZ20 = filter_CONDITION_S_deterministic(df_CLZ20)
        else:
            df_CLZ20 = filter_CONDITION_S(df_CLZ20)
    else:
        print("WARNING: skip_filter_S=True — the S-nonempty filter has been "
              "skipped. This is for testing only; output may include curves "
              "that should be filtered out and is not suitable for production.")

    df_CLZ20['source'] = 'CLZ20'

    # --------------------------------------------------------------
    # using criteria from [Zha16, Theorem 1.1 - 1.9]
    # --------------------------------------------------------------
    df_Zha16 = df
    df_Zha16 = merge_mwbsd(df_Zha16)
    # CLZ convention: L^(alg)(E,1) = L(E,1) / Omega, where Omega is the LMFDB
    # real_period. Under this convention Zhai's piecewise check (ord_2(L/Omega+)
    # = 0 for Delta<0, = 1 for Delta>0) collapses to ord_2(L^(alg)) = 0 in both
    # cases, since Omega = #real_components * Omega+ contributes the missing
    # factor of 2 to the 2-adic valuation when Delta>0.
    df_Zha16['L_alg'] = df_Zha16['special_value'] / df_Zha16['real_period']
    df_Zha16['L_alg_ord_2'] = df_Zha16['L_alg'].apply(lambda x: QQ(RR(x)).valuation(2))

    # E(Q)[2] = 0 & ord_2(L^(alg)) = 0 (Zhai Thm 1.1-1.9, CLZ convention)
    df_Zha16 = df_Zha16[~df_Zha16['torsion_primes'].apply(lambda x: 2 in x)]
    df_Zha16 = df_Zha16[df_Zha16['L_alg_ord_2'] == 0]
    df_Zha16['source'] = 'Zha16_no_2_tors'  # watermark
    # --------------------------------------------------------------

    # combine the results from CLZ20 and Zha16
    df = pd.concat([df_CLZ20, df_Zha16]).sort_values(by=['source','conductor']).reset_index(drop=True)
    df = df[['source','conductor','lmfdb_label','Clabel','ainvs','torsion_structure','sha']
            ].drop_duplicates(subset=['lmfdb_label'])
    res = list(zip(df['Clabel'].tolist(), df['source'].tolist(),
                   df['lmfdb_label'].tolist(), df['ainvs'].tolist(),
                   df['torsion_structure'].tolist(), df['sha'].tolist()))

    # Final gate: unconditionally check the 2-part of BSD on each surviving curve.
    # Curves that pass go to ec_file; curves where every 2-descent backend fails to
    # pin down ord_2(#Sha) go to fail_file for follow-up.
    # LMFDB curves are already optimal and torsion_structure / sha are already in
    # hand, so check_BSD_at_2 takes the precomputed two-torsion rank and ord_2(sha_an).
    passed, failed = [], []
    n_total = len(res)
    if skip_BSD_at_2_check:
        print(
            "WARNING: skip_BSD_at_2_check=True — the BSD(E,2) verification "
            "(check_BSD_at_2) is being skipped. Output curves are NOT verified "
            f"for BSD(E,2); they only satisfy the upstream filters ({n_total} "
            "curves). Use this only when BSD(E,2) for these curves has already "
            "been verified in a prior run.",
            flush=True,
        )
        passed = [(c, s, l) for c, s, l, _, _, _ in res]
    else:
        print(f"Running check_BSD_at_2 on {n_total} curves...", flush=True)
        for i, (cremona_label, source, lmfdb_label, ainvs, torsion_structure, sha_an) in enumerate(res, 1):
            E = EllipticCurve(ainvs)
            two_tor_rk = two_torsion_rank_from_structure(torsion_structure)
            sha_an_ord_2 = ZZ(sha_an).valuation(2)
            if check_BSD_at_2(E, two_tor_rk, sha_an_ord_2):
                passed.append((cremona_label, source, lmfdb_label))
            else:
                failed.append((cremona_label, source, lmfdb_label))
            if i % 100 == 0 or i == n_total:
                print(f"  {i}/{n_total} ({len(passed)} pass, {len(failed)} fail)", flush=True)

    # Derive a sibling fail file from EC_FILE: foo.txt -> foo_fail.txt
    base, ext = os.path.splitext(ec_file)
    fail_file = f"{base}_fail{ext or '.txt'}"

    # calculate the run time
    current_time = get_current_time_str()
    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes = int(elapsed_time // 60)
    seconds = int(elapsed_time % 60)

    bsd_check_note = (
        "BSD(E,2) check: SKIPPED (output not verified for BSD(E,2))"
        if skip_BSD_at_2_check
        else "BSD(E,2) check: run (output verified)"
    )
    header = (
        f"{current_time}\n"
        f"{run_summary}\n"
        f"{bsd_check_note}\n"
        f"Run took: {minutes} minutes {seconds} seconds\n\n"
        "cremona_label,source,lmfdb_label\n"
    )

    with open(ec_file, 'w') as f:
        f.write(header)
        for cremona_label, source, lmfdb_label in passed:
            f.write(f"{cremona_label}, {source}, {lmfdb_label}\n")

    if skip_BSD_at_2_check:
        print(f"SUCCESS!!! {len(passed)} curves output (BSD(E,2) check skipped), saved to {ec_file}.")
    elif failed:
        with open(fail_file, 'w') as f:
            f.write(header)
            for cremona_label, source, lmfdb_label in failed:
                f.write(f"{cremona_label}, {source}, {lmfdb_label}\n")
        print(f"SUCCESS!!! {len(passed)} curves passed BSD(E,2) check, saved to {ec_file}.")
        print(f"{len(failed)} curves failed, saved to {fail_file}.")
    else:
        # No failures this run; remove any stale fail file so it doesn't misrepresent state.
        if os.path.exists(fail_file):
            os.remove(fail_file)
        print(f"SUCCESS!!! {len(passed)} curves passed BSD(E,2) check, saved to {ec_file}. No failures.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Find elliptic curves with infinite families of quadratic twists satisfying BSD."
    )
    parser.add_argument(
        "--cond_upper_bound",
        type=int,
        default=None,
        help="Upper bound for conductor (e.g., 150 for testing)"
    )
    parser.add_argument(
        "--cond_lower_bound",
        type=int,
        default=None,
        help="Lower bound for conductor"
    )
    parser.add_argument(
        "--skip_filter_S",
        action="store_true",
        help="FOR TESTING ONLY. Skip the filter_CONDITION_S filter on the CLZ20 "
             "branch. This is unsound: skipping S lets through curves that should "
             "be filtered out. Off by default; do not enable for production runs."
    )
    parser.add_argument(
        "--use_probabilistic_S",
        action="store_true",
        help="Use the probabilistic S filter (filter_CONDITION_S) instead of the "
             "deterministic one (filter_CONDITION_S_deterministic, the default and "
             "the form stated in the paper)."
    )
    parser.add_argument(
        "--ec_file",
        type=str,
        default="output/ec_labels.txt",
        help="Path to save the output file"
    )
    parser.add_argument(
        "--skip_BSD_at_2_check",
        action="store_true",
        help="Skip the unconditional BSD(E,2) verification step (check_BSD_at_2). "
             "Default behavior is to run the check. Use this flag for faster re-runs "
             "when BSD(E,2) has already been verified on the same curves; output curves "
             "are then guaranteed only to satisfy the upstream filters, NOT BSD(E,2)."
    )
    args = parser.parse_args()

    bsd_infinite_twists(
        cond_upper_bound=args.cond_upper_bound,
        cond_lower_bound=args.cond_lower_bound,
        EC_FILE=args.ec_file,
        skip_filter_S=args.skip_filter_S,
        use_deterministic_S=(not args.use_probabilistic_S),
        skip_BSD_at_2_check=args.skip_BSD_at_2_check,
    )