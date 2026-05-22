"""Algorithm2.py

Computes admissible BSD twists for elliptic curves identified by Algorithm1.

This script reads elliptic curve labels from the output of Algorithm 1 and computes
the admissible quadratic twists M (up to a bound B) for which the full BSD conjecture
formula can be verified.

Two methods are implemented based on the source paper:
- CLZ20: Cai-Li-Zhai (2020) criteria for curves with E(Q)[2] = Z/2Z
- Zha16: Zhai (2016) criteria for curves without 2-torsion

Usage:
    sage -python Algorithm2.py <input_file>

    where <input_file> is a text file containing elliptic curve labels from Algorithm 1
    (e.g., ec_labels_150.txt).

    The output will be written to output/twists_of_<input_filename>.json
    For example:
        sage -python Algorithm2.py output/ec_labels_150.txt
    will produce:
        output/twists_of_ec_labels_150.json

References:
    [CLZ20] L. Cai, C. Li and S. Zhai, "On the 2-part of the Birch and
            Swinnerton-Dyer conjecture for quadratic twists of elliptic curves",
            J. Lond. Math. Soc. (2) 101 (2020), no. 2, 714-734.
    [Zha16] S. Zhai, "Non-vanishing theorems for quadratic twists of elliptic curves",
            Asian J. Math. 20 (2016), no. 3, 475-502.
"""

import sys
import os
# Add repo root to path (go up 2 levels: infinite_bsd -> ants_xvii -> repo_root)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import json
from sage.all import (
    EllipticCurve, kronecker_symbol, ZZ, gcd,
    prime_range, fundamental_discriminant, NumberField
)
from lmfdb import db

# =============================================================================
# CONFIGURATION
# =============================================================================

# Bound for admissible twists M
TWIST_BOUND = 1000

# LMFDB database handle
ecq = db.ec_curvedata

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def p_inert_in_F(p: int, F: NumberField) -> bool:
    """
    Check if a prime p is inert in the number field F.

    A prime is inert if it remains prime in the ring of integers of F,
    i.e., the ideal (p) is a prime ideal in O_F.

    Args:
        p: A prime number
        F: A number field

    Returns:
        True if p is inert in F, False otherwise
    """
    I = F.ideal(p)
    return I.is_prime()


def disc_valuation_condition(primes_dividing_M: list,
                             conductor_primes: list,
                             disc_valuations: dict) -> bool:
    """
    Test the ramification-style condition

        for every p | M, there exists q | N, q != p,
        such that p does not divide ord_q(Delta_E).

    Mirrors `is_ramified` in Algorithm1.py but with the outer loop ranging
    over primes of the twist M rather than primes of the base conductor N.

    Args:
        primes_dividing_M: primes p with p | M.
        conductor_primes: primes q with q | N (conductor of E).
        disc_valuations: dict mapping each q in conductor_primes
            to ord_q(Delta_E).

    Returns:
        True iff the condition holds.
    """
    return all(
        any(disc_valuations[q] % p != 0 for q in conductor_primes if q != p)
        for p in primes_dividing_M
    )


# =============================================================================
# ADMISSIBLE TWIST FUNCTIONS
# =============================================================================

def get_admissible_twists_CLZ(E: EllipticCurve, B: int = 150) -> list:
    """
    Compute admissible BSD twists using the Cai-Li-Zhai (2020) criteria.

    For curves E with E(Q)[2] = Z/2Z, this finds squarefree integers M such that
    the quadratic twist E_M satisfies the conditions of [CLZ20, Theorem 1.5].

    Conditions on M:
        (a) M is squarefree and gcd(M, N) = 1 where N is the conductor
        (b) a_p(E) is not divisible by p for all primes p | M
        (c) p ≡ 1 (mod 4) and ord_2(a_p) = 1 for all primes p | M
        (d) M ≡ 1 (mod 8) and (disc(M)/q) = 1 for all odd primes q | N

    Args:
        E: An elliptic curve over Q
        B: Upper bound for twists to consider (default: 150)

    Returns:
        List of admissible twist values M
    """
    conductor = E.conductor()
    conductor_primes = ZZ(conductor).prime_divisors()

    # Precompute a_p values for efficiency
    a_p_dict = {p: E.ap(p) for p in prime_range(B)}

    # Precompute ord_q(Delta_E) at each bad prime q (used by the
    # disc-valuation ramification condition below)
    disc_valuations = {q: E.discriminant().valuation(q) for q in conductor_primes}

    admissible_twists = []

    for M in range(1, B + 1):

        # Condition (a): M is squarefree
        if not ZZ(M).is_squarefree():
            continue

        # Condition (a): gcd(M, conductor) = 1
        if gcd(M, conductor) != 1:
            continue

        primes_dividing_M = ZZ(M).prime_divisors()

        # Ramification condition: for every p | M, there exists q | N, q != p,
        # such that p does not divide ord_q(Delta_E).
        if not disc_valuation_condition(primes_dividing_M, conductor_primes, disc_valuations):
            continue

        # Condition (b): a_p not divisible by p for all p | M
        if not all(a_p_dict[p] % p != 0 for p in primes_dividing_M):
            continue

        # Condition (c): p ≡ 1 (mod 4) and ord_2(a_p) = 1 for all p | M
        if not all((p % 4 == 1) and ((p+1-a_p_dict[p]).valuation(2) == 1) for p in primes_dividing_M):
            continue

        # Condition (d): M ≡ 1 (mod 8) and Kronecker conditions
        if M % 8 != 1:
            continue

        disc_M = fundamental_discriminant(M)
        if all(kronecker_symbol(disc_M, p) == 1 for p in conductor_primes if p != 2):
            admissible_twists.append(M)

    return admissible_twists


def get_admissible_twists_Zhai(E: EllipticCurve, B: int = 150) -> list:
    """
    Compute admissible BSD twists using the Zhai (2016) criteria.

    For curves E without 2-torsion, this finds squarefree integers M such that
    the quadratic twist E_M satisfies the conditions of [Zha16, Theorems 1.1-1.9].

    Conditions on M:
        (a) M is squarefree and gcd(M, N) = 1 where N is the conductor
        (b) a_p(E) is not divisible by p for all primes p | M
        (c) M ≡ 1 (mod 4)
        (d) All primes p | M are inert in Q(E[2]) (the 2-division field)
        (e) Kronecker symbol conditions on conductor primes

    Note: If disc(E) > 0, only positive M are considered.

    Args:
        E: An elliptic curve over Q
        B: Upper bound for twists to consider (default: 150)

    Returns:
        List of admissible twist values M
    """
    conductor = E.conductor()
    conductor_primes = ZZ(conductor).prime_divisors()

    # Precompute a_p values for efficiency
    a_p_dict = {p: E.ap(p) for p in prime_range(B)}

    # Precompute ord_q(Delta_E) at each bad prime q (used by the
    # disc-valuation ramification condition below)
    disc_valuations = {q: E.discriminant().valuation(q) for q in conductor_primes}

    admissible_twists = []

    # Precompute the 2-division field (only need to do this once per curve)
    two_division_poly = E.two_division_polynomial()
    two_division_field = NumberField(two_division_poly, 'a')

    for M in range(-B, B + 1):
        if M == 0:
            continue

        # If discriminant is positive, only consider positive M
        if E.discriminant() > 0 and M < 0:
            continue

        # Condition (a): M is squarefree
        if not ZZ(M).is_squarefree():
            continue

        # Condition (a): gcd(M, conductor) = 1
        if gcd(M, conductor) != 1:
            continue

        primes_dividing_M = ZZ(M).prime_divisors()

        # Ramification condition: for every p | M, there exists q | N, q != p,
        # such that p does not divide ord_q(Delta_E).
        if not disc_valuation_condition(primes_dividing_M, conductor_primes, disc_valuations):
            continue

        # Condition (b): a_p not divisible by p for all p | M
        if not all(a_p_dict[p] % p != 0 for p in primes_dividing_M):
            continue

        # Condition (c): M ≡ 1 (mod 4)
        if M % 4 != 1:
            continue

        # Condition (d): all primes dividing M are inert in Q(E[2])
        if not all(p_inert_in_F(p, two_division_field) for p in primes_dividing_M):
            continue

        # Condition (e): Kronecker symbol conditions
        disc_M = fundamental_discriminant(M)
        if all(kronecker_symbol(disc_M, p) == 1 for p in conductor_primes):
            admissible_twists.append(M)

    return admissible_twists


# =============================================================================
# MAIN
# =============================================================================

def main():
    """
    Main entry point: read curve labels, compute admissible twists, save results.
    """
    if len(sys.argv) != 2:
        print("Usage: sage -python Algorithm2.py <input_file>")
        print("Example: sage -python Algorithm2.py output/ec_labels_150.txt")
        sys.exit(1)

    input_file = sys.argv[1]

    # Generate output filename from input filename
    input_basename = os.path.basename(input_file)
    input_name_without_ext = os.path.splitext(input_basename)[0]
    output_file = f"output/twists_of_{input_name_without_ext}.json"

    print(f"Computing admissible twists up to bound B = {TWIST_BOUND}")
    print(f"Reading labels from: {input_file}")

    # Read elliptic curve labels from Algorithm 1 output
    with open(input_file, 'r') as file:
        labels = [line.strip() for line in file.readlines()]

    # Algorithm1 writes metadata, a blank line, then a "cremona_label,source,lmfdb_label"
    # column header, then data. Accept only lines that look like a data row: at
    # least 3 comma-separated parts where the second part is a known source tag.
    KNOWN_SOURCES = {'CLZ20', 'Zha16_no_2_tors'}
    rows = []
    for label_line in labels:
        parts = [p.strip() for p in label_line.split(',')]
        if len(parts) < 3 or parts[1] not in KNOWN_SOURCES:
            continue
        rows.append((parts[0], parts[1], parts[2]))  # (cremona, source, lmfdb)

    # Bulk-fetch ainvs in chunks instead of one ecq.lookup per curve.
    # For ~40k curves this turns ~40k DB round-trips into ~40 batched ones.
    all_labels = [r[2] for r in rows]
    print(f"Bulk-fetching ainvs for {len(all_labels)} curves from LMFDB...")
    ainvs_by_label = {}
    BATCH = 1000
    for start in range(0, len(all_labels), BATCH):
        chunk = all_labels[start:start + BATCH]
        for r in ecq.search({'lmfdb_label': {'$in': chunk}},
                            projection=['lmfdb_label', 'ainvs']):
            ainvs_by_label[r['lmfdb_label']] = r['ainvs']
    print(f"  retrieved ainvs for {len(ainvs_by_label)} curves")

    results = {}
    for cremona_label, source, lmfdb_label in rows:
        ainvs = ainvs_by_label.get(lmfdb_label)
        if ainvs is None:
            print(f"  WARNING: ainvs not found for {lmfdb_label}; skipping")
            continue
        E = EllipticCurve(ainvs)

        # Compute admissible twists based on source paper
        if source == 'Zha16_no_2_tors':
            twists = get_admissible_twists_Zhai(E, TWIST_BOUND)
        elif source == 'CLZ20':
            twists = get_admissible_twists_CLZ(E, TWIST_BOUND)
        else:
            raise NotImplementedError(f"Source {source} not implemented yet.")

        results[cremona_label] = twists
        print(f"  {cremona_label} ({source}): {len(twists)} admissible twists")

    # Save results to JSON
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=4)

    print(f"\nResults saved to: {output_file}")


if __name__ == '__main__':
    main()
