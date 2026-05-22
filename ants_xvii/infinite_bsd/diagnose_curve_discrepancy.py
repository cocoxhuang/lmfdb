"""diagnose_curve_discrepancy.py

For four Cremona-labeled curves (62a1, 66b1, 105a1, 141c1) — which
Zhai's paper claims pass his theorem's conditions but our Algorithm1.py
rejects at conductor < 150 — run every filter from Algorithm1's pipeline
independently and report PASS / FAIL / SKIP / INFO per condition, with
the specific value when relevant.

Output: per-curve report to stdout and to output/discrepancy_report.txt.

Usage:
    sage -python diagnose_curve_discrepancy.py
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import pandas as pd
from lmfdb import db
from sage.all import EllipticCurve, Integer, QQ, RR, ZZ, PolynomialRing

from Algorithm1 import (
    is_ramified,
    check_BSD_at_2,
    two_torsion_rank_from_structure,
    filter_CONDITION_S,
)

# bound used by Algorithm1's probabilistic S filter
S_BOUND = 10000

CURVES = ['62a1', '66b1', '105a1', '141c1']

ecq = db.ec_curvedata
ec_classdata = db.ec_classdata
ec_mwbsd = db.ec_mwbsd

ECQ_COLS = ['ainvs', 'lmfdb_label', 'Clabel', 'conductor', 'rank',
            'absD', 'bad_primes', 'isogeny_degrees', 'semistable',
            'num_bad_primes', 'optimality', 'torsion_structure',
            'sha', 'lmfdb_iso']
MWBSD_COLS = ['lmfdb_label', 'real_period', 'special_value', 'tamagawa_product']


def fetch_data():
    """One bulk query per LMFDB table for all four curves."""
    rows = list(ecq.search({'Clabel': {'$in': CURVES}}, projection=ECQ_COLS))
    by_clabel = {r['Clabel']: r for r in rows}

    iso_labels = list({r['lmfdb_iso'] for r in rows})
    cd_rows = list(ec_classdata.search({'lmfdb_iso': {'$in': iso_labels}},
                                       projection=['lmfdb_iso', 'aplist']))
    a3_by_iso = {r['lmfdb_iso']: r['aplist'][1] for r in cd_rows}

    lmfdb_labels = [r['lmfdb_label'] for r in rows]
    bsd_rows = list(ec_mwbsd.search({'lmfdb_label': {'$in': lmfdb_labels}},
                                    projection=MWBSD_COLS))
    bsd_by_lmfdb = {r['lmfdb_label']: r for r in bsd_rows}

    return by_clabel, a3_by_iso, bsd_by_lmfdb


def emit(out, line):
    out.append(line)
    print(line)


def status(passed):
    return "PASS" if passed else "FAIL"


def diagnose_curve(clabel, row, a3, bsd_row, out):
    ainvs = row['ainvs']
    lmfdb_label = row['lmfdb_label']
    N = row['conductor']
    bad_primes = list(row['bad_primes'])
    isogeny_degrees = list(row['isogeny_degrees'])
    torsion_structure = row['torsion_structure']
    minimal_disc = Integer(row['absD'])

    E = EllipticCurve(ainvs)
    delta_E = E.discriminant()

    failures = []

    emit(out, "")
    factorization = '·'.join(str(p) for p in bad_primes)
    emit(out, f"=== {clabel}  (LMFDB {lmfdb_label})  ainvs={ainvs}  "
              f"N={N}={factorization}  rank={row['rank']} ===")
    emit(out, "")

    # 1. Squarefree, not prime
    ok = bool(row['semistable']) and row['num_bad_primes'] >= 2
    emit(out, f"  [{status(ok)}] N squarefree, not prime          "
              f"semistable={row['semistable']}, num_bad_primes={row['num_bad_primes']}")
    if not ok:
        failures.append("N not squarefree or N prime")

    # 2. Isogeny condition at A
    A = set()
    if 3 in bad_primes or abs(a3) == 3:
        A.add(3)
    if 5 in bad_primes:
        A.add(5)
    if 7 in bad_primes:
        A.add(7)
    isog_primes = {p for p in isogeny_degrees if Integer(p).is_prime()}
    overlap = A & isog_primes
    ok = not overlap
    emit(out, f"  [{status(ok)}] No p-isogeny for p in A          "
              f"A={sorted(A)}, prime-degree isogenies={sorted(isog_primes)}, a3={a3}")
    if not ok:
        failures.append(f"rational p-isogeny exists for p ∈ A ∩ isog_primes = {sorted(overlap)}")

    # 3. Ramification
    ok = is_ramified(bad_primes, minimal_disc)
    witnesses = []
    for p in bad_primes:
        found = None
        for q in bad_primes:
            if q != p and minimal_disc.valuation(q) % p != 0:
                found = (q, int(minimal_disc.valuation(q)))
                break
        witnesses.append((p, found))
    detail = '; '.join(
        f"p={p}→q={w[0]},ord_q(Δ)={w[1]}" if w else f"p={p}: NO WITNESS"
        for p, w in witnesses
    )
    emit(out, f"  [{status(ok)}] Ramification                    {detail}")
    if not ok:
        failures.append("ramification condition fails for some prime of N")

    # 4. Optimal
    ok = row['optimality'] == 1
    emit(out, f"  [{status(ok)}] Optimal                          optimality={row['optimality']}")
    if not ok:
        failures.append("not optimal")

    # 5. rank 0
    ok = row['rank'] == 0
    emit(out, f"  [{status(ok)}] rank(E)=0                        rank={row['rank']}")
    if not ok:
        failures.append(f"rank {row['rank']} ≠ 0")

    # 6. 2-torsion branch
    two_tor_rk = two_torsion_rank_from_structure(torsion_structure)
    if two_tor_rk == 1:
        branch = "CLZ20"
    elif two_tor_rk == 0:
        branch = "Zhai"
    else:
        branch = "NEITHER"
    emit(out, f"  [INFO] 2-torsion branch                  torsion_structure={torsion_structure}, "
              f"dim_F2 E(Q)[2]={two_tor_rk} → {branch} branch")
    if branch == "NEITHER":
        failures.append(f"E(Q)[2] dim={two_tor_rk}; pipeline only handles dim 0 or 1")

    # L^(alg) = special_value / real_period in CLZ convention
    real_period = bsd_row['real_period']
    special_value = bsd_row['special_value']
    L_alg = QQ(RR(special_value) / RR(real_period))
    L_alg_ord_2 = L_alg.valuation(2)

    # 7. CLZ20 branch checks
    if branch == "CLZ20":
        # 7a
        ok = L_alg_ord_2 == -1
        emit(out, f"  [{status(ok)}] ord_2(L^(alg)) = −1 (CLZ20)     "
                  f"actual ord_2={L_alg_ord_2}, L^(alg)={L_alg} ≈ {float(L_alg):.6g}")
        if not ok:
            failures.append(f"CLZ20: ord_2(L^(alg))={L_alg_ord_2}, expected −1")

        # 7b/c: build E' from the (unique) order-2 generator of E(Q)_tors
        E_gens = E.torsion_subgroup().gens()
        if len(E_gens) != 1:
            emit(out, f"  [SKIP] E' analysis              E(Q)_tors has {len(E_gens)} generators (not cyclic)")
        else:
            P = E_gens[0]
            order = P.order()
            if order % 2 != 0:
                emit(out, f"  [SKIP] E' analysis              E(Q)_tors order {order} odd; no E(Q)[2]")
            else:
                E_two_gen = (order // 2) * P
                C = E(E_two_gen)
                E_prime = E.isogeny_codomain(C)
                Ep_invariants = E_prime.torsion_subgroup().invariants()
                Ep_tors_order = E_prime.torsion_order()
                n_gen = E_prime.torsion_subgroup().ngens()

                # 7b
                ok = (n_gen == 1 and Ep_tors_order % 2 == 0)
                emit(out, f"  [{status(ok)}] E'(Q)[2] ≅ Z/2Z (CLZ20)         "
                          f"E'(Q)_tors invariants={Ep_invariants}, order={Ep_tors_order}")
                if not ok:
                    failures.append(f"E'(Q)[2] ≠ Z/2Z (E'(Q)_tors invariants={Ep_invariants})")

                # 7c — Sha(E')[2] = 0 via descent
                try:
                    Ep_sha_an = E_prime.sha().an().round()
                except Exception as e:
                    try:
                        Ep_sha_an = ecq.lucky(
                            {'ainvs': [int(x) for x in E_prime.ainvs()]},
                            projection='sha')
                    except Exception:
                        Ep_sha_an = None
                if Ep_sha_an is None:
                    emit(out, f"  [ERR ] Sha(E')[2] = 0 (CLZ20)         could not obtain sha_an(E')")
                    failures.append("could not obtain sha_an(E')")
                else:
                    sha_ord_2_Ep = ZZ(Ep_sha_an).valuation(2)
                    ok = check_BSD_at_2(E_prime, two_tor_rk=1, sha_an_ord_2=sha_ord_2_Ep)
                    emit(out, f"  [{status(ok)}] Sha(E')[2] = 0 (CLZ20)           "
                              f"sha_an(E')={Ep_sha_an} (ord_2={sha_ord_2_Ep}), descent verdict={ok}")
                    if not ok:
                        failures.append(f"Sha(E')[2] = 0 not verified by descent (sha_an(E')={Ep_sha_an})")

        # 7d — square check on f'(x_0), -f'(x_0), -Δ_E (deterministic S filter)
        try:
            E_sw = E.short_weierstrass_model()
            A_coef = E_sw.a4()
            B_coef = E_sw.a6()
            P_ring = PolynomialRing(QQ, 'x')
            f = P_ring([B_coef, A_coef, 0, 1])
            roots = f.roots()
            if len(roots) != 1:
                emit(out, f"  [SKIP] square check (CLZ20)      f(x) has {len(roots)} rational roots, expected exactly 1")
            else:
                x0 = roots[0][0]
                fp_x0 = 3 * x0**2 + A_coef
                # in short Weierstrass: -Δ_E = 16(4A^3 + 27B^2)
                neg_delta_sw = 4 * A_coef**3 + 27 * B_coef**2  # equal to -Δ_E / 16
                checks = [
                    ("f'(x_0)",      fp_x0,            fp_x0.is_square()),
                    ("-f'(x_0)",     -fp_x0,           (-fp_x0).is_square()),
                    ("-Δ_E/16",      neg_delta_sw,     neg_delta_sw.is_square()),
                ]
                any_square = any(c[2] for c in checks)
                ok = not any_square
                detail = ', '.join(
                    f"{name}={val} ({'square' if sq else 'not square'})"
                    for name, val, sq in checks
                )
                emit(out, f"  [{status(ok)}] no square in {{f'(x_0),−f'(x_0),−Δ_E}} (CLZ20)  {detail}")
                if not ok:
                    squared = [name for name, _, sq in checks if sq]
                    failures.append(f"square check failed: {squared} is/are squares in Q")
        except Exception as e:
            emit(out, f"  [ERR ] square check (CLZ20)      {e}")

        # 7e — probabilistic S-nonempty filter (independent cross-check of 7d).
        # S = { q ≡ 1 (mod 4) : q ∤ N, ord_2(#E(F_q)) = 1 }; need at least one
        # witness q ≤ S_BOUND. Calls filter_CONDITION_S directly on a single-row
        # DataFrame so the diagnostic uses the exact same code path as the pipeline.
        df_single = pd.DataFrame([{'ainvs': ainvs, 'conductor': N}])
        S_survivors = filter_CONDITION_S(df_single, S_bound=S_BOUND)
        # find an explicit witness if one exists, for the report
        from sage.all import primes
        witness = None
        for q in primes(S_BOUND):
            if q % 4 == 1 and N % q != 0:
                N_q = ZZ(1 + q - E.ap(q))
                if N_q.valuation(2) == 1:
                    witness = q
                    break
        ok = len(S_survivors) > 0
        if ok:
            detail = f"witness q={witness}, ord_2(#E(F_q))=1"
        else:
            detail = f"no q ≡ 1 mod 4 with q ∤ N and ord_2(#E(F_q))=1 found up to {S_BOUND}"
        emit(out, f"  [{status(ok)}] S non-empty (CLZ20, probabilistic) {detail}")
        if not ok:
            failures.append(f"S empty up to {S_BOUND} (probabilistic check confirms deterministic)")

    elif branch == "Zhai":
        ok = L_alg_ord_2 == 0
        emit(out, f"  [{status(ok)}] ord_2(L^(alg)) = 0 (Zhai)        "
                  f"actual ord_2={L_alg_ord_2}, L^(alg)={L_alg} ≈ {float(L_alg):.6g}")
        if not ok:
            failures.append(f"Zhai: ord_2(L^(alg))={L_alg_ord_2}, expected 0")

    # 9. Final check_BSD_at_2 on E itself
    sha = row['sha']
    sha_ord_2 = ZZ(sha).valuation(2)
    ok = check_BSD_at_2(E, two_tor_rk=two_tor_rk, sha_an_ord_2=sha_ord_2)
    emit(out, f"  [{status(ok)}] check_BSD_at_2(E)                "
              f"sha_an(E)={sha} (ord_2={sha_ord_2}), descent verdict={ok}")
    if not ok:
        failures.append(f"check_BSD_at_2(E) failed (sha_an(E)={sha})")

    emit(out, "")
    if failures:
        emit(out, f"  ⇒ {clabel} is REJECTED by Algorithm1 because:")
        for f in failures:
            emit(out, f"        • {f}")
    else:
        emit(out, f"  ⇒ {clabel} passes ALL conditions — should appear in Algorithm1 output.")


def main():
    print(f"Diagnosing {len(CURVES)} curves: {CURVES}")
    print("Fetching LMFDB data (3 bulk queries)...")
    by_clabel, a3_by_iso, bsd_by_lmfdb = fetch_data()

    out = []
    emit(out, "Diagnostic report: Zhai-list discrepancy")
    emit(out, "=" * 80)
    emit(out, f"Curves: {CURVES}")

    for clabel in CURVES:
        if clabel not in by_clabel:
            emit(out, "")
            emit(out, f"!!! {clabel}: not found in LMFDB ec_curvedata.")
            continue
        row = by_clabel[clabel]
        a3 = a3_by_iso[row['lmfdb_iso']]
        bsd_row = bsd_by_lmfdb[row['lmfdb_label']]
        diagnose_curve(clabel, row, a3, bsd_row, out)

    report_path = "output/discrepancy_report.txt"
    with open(report_path, 'w') as f:
        f.write('\n'.join(out) + '\n')
    print()
    print(f"Report saved to {report_path}")


if __name__ == "__main__":
    main()
