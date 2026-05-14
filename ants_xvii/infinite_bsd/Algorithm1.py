"""Algorithm1.py

This script generates a text file which includes elliptic curves with
a family of quadratic twists that satisfy the full BSD conjecture formula.
It corresponds to Algorithm 1 in the paper.

To run this, execute the following command from this level of the directory in the terminal:

    sage -python Algorithm1.py

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
                    ZZ, primes, fundamental_discriminant, gcd, kronecker_symbol)  
import time
from datetime import datetime

# Load the elliptic curve table from the LMFDB
ecq = db.ec_curvedata
ec_classdata = db.ec_classdata
ec_mwbsd = db.ec_mwbsd

# Information columns from lmfdb
ECQ_COLS = ['ainvs', 'lmfdb_label', 'conductor', 'rank', 'torsion', 'absD', 'bad_primes', 'manin_constant', 
'regulator', 'sha', 'lmfdb_iso', 'torsion_structure', 'torsion_primes', 'signD', 'Clabel']
CLASSDATA_COLS = ['lmfdb_iso', 'aplist']
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

def ordinary_357(a_ps: list[int]) -> list[int]:
    '''returns the ordinary primes among 3,5,7 for a given elliptic curve'''
    primes = [3,5,7]
    return [primes[i] for i in range(len(primes)) if a_ps[i] % primes[i] != 0]

def get_a3_a5_a7(df: pd.DataFrame) -> pd.DataFrame:
    '''returns a dataframe with a3,a5,a7 columns merged from classdata'''
    classdata_query = {'lmfdb_iso': {'$in': df['lmfdb_iso'].tolist()}}
    classdata_payload = ec_classdata.search(classdata_query, projection=CLASSDATA_COLS)
    classdata_df = pd.DataFrame(list(classdata_payload))
    classdata_df['a3'] = classdata_df['aplist'].apply(lambda x: x[1])
    classdata_df['a5'] = classdata_df['aplist'].apply(lambda x: x[2])
    classdata_df['a7'] = classdata_df['aplist'].apply(lambda x: x[3])
    classdata_df.drop(columns=['aplist'], inplace=True)
    return classdata_df

def merge_mwbsd(df: pd.DataFrame) -> pd.DataFrame:
    '''merges mwbsd data into the given dataframe'''
    mwbsd_query = {'lmfdb_label': {'$in': df['lmfdb_label'].tolist()}}
    mwbsd_payload = ec_mwbsd.search(mwbsd_query, projection=MWBSD_COLS)
    mwbsd_df = pd.DataFrame(list(mwbsd_payload))
    df = pd.merge(df, mwbsd_df, on='lmfdb_label', how='inner')
    return df

def filter_CONDITION_p_isogeny(df: pd.DataFrame) -> pd.DataFrame:
    '''CONDITION: no rational odd prime isogenies'''
    data_idx = []
    for index,curve in df.iterrows():
        ainvs = curve['ainvs']
        bad_primes = curve['bad_primes']

        # CONDITION: irreducibility CONDITION at certain odd primes
        CONDITION = True
        a_3, a_5, a_7 = curve['a3'], curve['a5'], curve['a7']
        non_isogeny_primes = set(bad_primes) and set([3,5,7])
        # non_isogeny_primes.update(ordinary_357([a_3, a_5, a_7]))
    
        E = EllipticCurve(ainvs)
        isogeny_primes = [phi.degree() for phi in E.isogenies_prime_degree()]
        CONDITION = (all(p not in isogeny_primes for p in non_isogeny_primes))

        if CONDITION:
            data_idx.append(index)

    return df.loc[data_idx]

def filter_CONDITION_no_3_isogeny(df: pd.DataFrame) -> pd.DataFrame:
    '''CONDITION: no rational 3 isogenies'''
    data_idx = []
    for index,curve in df.iterrows():
        ainvs = curve['ainvs']

        # CONDITION: 3 is not a isogeny prime
        E = EllipticCurve(ainvs)
        isogeny_primes = [phi.degree() for phi in E.isogenies_prime_degree()]
        CONDITION = 3 not in isogeny_primes

        if CONDITION:
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
    '''CONDITION: sha(E')[2] = 1 and E'(Q)[2] = Z/2Z, where E' is the 2-isogenous curve of E'''
    data_idx = []
    for index,curve in df.iterrows():
        ainvs = curve['ainvs']
        E = EllipticCurve(ainvs)

        # find the generator of the 2-torsion subgroup
        E_torsion_gens = E.torsion_subgroup().gens()
        assert len(E_torsion_gens) == 1
        P = E_torsion_gens[0]
        if P.order() == 2:
            E_two_torsion_gen = P
        elif P.order() == 4:
            E_two_torsion_gen = 2 * P
        elif P.order() == 6:
            E_two_torsion_gen = 3 * P
        elif P.order() == 8:
            E_two_torsion_gen = 4 * P
        elif P.order() == 10:
            E_two_torsion_gen = 5 * P
        else:
            raise ValueError("Unexpected torsion order")
        assert E_two_torsion_gen.order() == 2   
        
        # Get sha(E')
        C = E(E_two_torsion_gen)
        E_prime = E.isogeny_codomain(C)
        try:
            E_prime_sha_order = E_prime.sha().an().round()  # the "round" is for the case where analytic rank is > 1
        except Exception as e:
            print(f"Error calculating sha for curve {curve['lmfdb_label']}: {e}")
            print("Getting sha value from db")
            try:
                E_prime_sha_order = ecq.lucky({'ainvs':[int(x) for x in E_prime.ainvs()]}, projection='sha')
            except:
                import pdb; pdb.set_trace()

        # if E_prime_sha_order == 1:
        # E'(Q)[2] = Z/2Z, on the premise that E(Q)[2] = Z/2Z
        if E_prime_sha_order % 2 != 0:
            n_gen = E_prime.torsion_subgroup().ngens()
            E_prime_tors_order = E_prime.torsion_order()
            if n_gen == 1 and E_prime_tors_order % 2 == 0:
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

def ord_2_two_torsion_order(torsion_structure, torsion):
    '''returns ord_2( the order of the 2-torsion subgroup)'''
    if len(torsion_structure) == 1:
        two_torsion_val_2 = min(ZZ(torsion).valuation(2), 1)
    else:
        two_torsion_val_2 = ZZ(torsion_structure[0]).valuation(2) + ZZ(torsion_structure[1]).valuation(2)
    return two_torsion_val_2

# Main function
def bsd_infinite_twists(cond_upper_bound:int = None, cond_lower_bound:int = None,
                        EC_FILE = 'output/ec_labels.txt',
                        full_2_torsion = False):
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
        full_2_torsion (bool): If True, considers curves with full 2-torsion instead of just Z/2Z torsion. 
            Defaults to False as there is currently no theoretical result for infinite BSD twists in this case.
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
                'manin_constant' : {'$mod': [1, 2]} # CONDITION 1f: manin constant is odd, using psycodict's ordering on '$mod
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

    # CONDITION: a3 in {-2, -1, 0, 1, 2}
    classdata_df = get_a3_a5_a7(df)
    lmfdb_iso_labels_a3_0_cond = classdata_df[classdata_df['a3'].isin({-2, -1, 0, 1, 2})]
    # df = pd.merge(df, lmfdb_iso_labels_a3_cond, on='lmfdb_iso', how='inner')

    # CONDITION: a3 in {-3, 3} and there is no rational 3-isogeny as an alternative to the above a3 condition
    # Reference: On the Iwasawa main conjectures for modular forms at non-ordinary primes (2018)
    lmfdb_iso_labels_alt_a3_cond = classdata_df[classdata_df['a3'].isin({-3, 3})]
    lmfdb_iso_labels_alt_a3_cond = pd.merge(lmfdb_iso_labels_alt_a3_cond, df[['lmfdb_iso', 'ainvs']], on='lmfdb_iso', how='inner')
    lmfdb_iso_labels_alt_a3_cond = filter_CONDITION_no_3_isogeny(lmfdb_iso_labels_alt_a3_cond)
    lmfdb_iso_labels_alt_a3_cond = lmfdb_iso_labels_alt_a3_cond[['lmfdb_iso']]

    lmfdb_iso_labels_a3_cond = pd.merge(lmfdb_iso_labels_a3_0_cond, lmfdb_iso_labels_alt_a3_cond, on='lmfdb_iso', how='outer')
    df = pd.merge(df, lmfdb_iso_labels_a3_cond, on='lmfdb_iso', how='inner')

    # CONDITION: no rational odd prime isogenies
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
    df_CLZ20 = filter_CONDITION_S(df_CLZ20)

    df_CLZ20['source'] = 'CLZ20'

    # --------------------------------------------------------------
    # using criteria from [Zha16, Theorem 1.1 - 1.9]
    # --------------------------------------------------------------
    df_Zha16 = df
    df_Zha16 = merge_mwbsd(df_Zha16)
    df_Zha16['real_components'] = df_Zha16['ainvs'].apply(lambda x: EllipticCurve(x).real_components())  
    df_Zha16['L_alg'] = (df_Zha16['special_value'] * df_Zha16['real_components']/(df_Zha16['real_period']))
    df_Zha16['L_alg_ord_2'] = df_Zha16['L_alg'].apply(lambda x: QQ(RR(x)).valuation(2))
    
    # part 1: NO 2 torsion & ( CONDITION on L_alg based on signD )
    df_Zha16_no_2_tors = df_Zha16[~df_Zha16['torsion_primes'].apply(lambda x: 2 in x)]   
    df_Zha16_no_2_tors = df_Zha16_no_2_tors.loc[( (df_Zha16_no_2_tors['L_alg_ord_2'] == 0) & (df_Zha16_no_2_tors['signD'] < 0) )| ( (df_Zha16_no_2_tors['L_alg_ord_2'] == 1) & (df_Zha16_no_2_tors['signD'] > 0) )]

    df_Zha16_no_2_tors['source'] = 'Zha16_no_2_tors' # watermark

    if not full_2_torsion:
        df_Zha16 = df_Zha16_no_2_tors
    else:   # part 2: 2 torsion & ( CONDITION on L_alg based on signD )
        df_Zha16_2_tors = df_Zha16[df_Zha16['torsion_primes'].apply(lambda x: 2 in x)]
        
        # L value conditions
        df_Zha16_2_tors = df_Zha16_2_tors.loc[( (df_Zha16_2_tors['L_alg_ord_2'] != 0) & (df_Zha16_2_tors['signD'] < 0) )| ( (df_Zha16_2_tors['L_alg_ord_2'] != 1) & (df_Zha16_2_tors['signD'] > 0) )]

        df_Zha16_2_tors['source'] = 'Zha16_2_tors' # watermark

        # combine the two
        df_Zha16 = pd.concat([df_Zha16_no_2_tors, df_Zha16_2_tors])
    # --------------------------------------------------------------

    # combine the results from CLZ20 and Zha16
    df = pd.concat([df_CLZ20, df_Zha16]).sort_values(by=['source','conductor']).reset_index(drop=True)
    df = df[['source','conductor','lmfdb_label','Clabel']].drop_duplicates(subset=['lmfdb_label'])
    lmfdb_labels = df['lmfdb_label'].tolist()
    cremona_labels = df['Clabel'].tolist()
    sources = df['source'].tolist()
    res = list(zip(cremona_labels, sources, lmfdb_labels))

    # calculate the run time
    current_time = get_current_time_str()
    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes = int(elapsed_time // 60)
    seconds = int(elapsed_time % 60)

    # save lmfdb labels
    with open(ec_file, 'w') as f:
        # Write the timestamp at the top of the file
        f.write(f"{current_time}\n")
        f.write(f"{run_summary}\n")
        f.write(f"Run took: {minutes} minutes {seconds} seconds\n\n")
        f.write("cremona_label,source,lmfdb_label\n")
        for cremona_label, source, lmfdb_label in res:
            f.write(f"{cremona_label}, {source}, {lmfdb_label}\n")


    print(f"SUCCESS!!! EC labels saved to {ec_file}.")

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
    args = parser.parse_args()

    bsd_infinite_twists(
        cond_upper_bound=args.cond_upper_bound,
        cond_lower_bound=args.cond_lower_bound
    )