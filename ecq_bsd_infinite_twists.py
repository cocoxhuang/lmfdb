"""ecq_bsd_infinite_twists.py

This script generates data files for elliptic curves with conductor less than 1024.
The data files contain the normalized a_p values for the first 50 primes,
as well as the conductor, rank, torsion, adelic level, adelic index, adelic genus,
and sha of the curve. The data files are saved in the parquet format.

To run this, execute the following command in the terminal:

    sage -python ecq_bsd_infinite_twists.py

"""


from lmfdb import db
import pandas as pd
import numpy as np

from sage.all import EllipticCurve, primes_first_n, round, Integer
import time

OUTPUT_FILE = 'data/output.txt'

# Load the elliptic curve table from the LMFDB
ecq = db.ec_curvedata
ec_localdata = db.ec_localdata
ec_mwbsd = db.ec_mwbsd

# Run a search query to get only the rank 0 or 1 curves.
# For now we set a limit to make things faster

NUM_AP_VALS = 100  # Number of primes to use for the a_p values
ADDITIONAL_COLS = ['conductor', 'rank', 'torsion', 'absD', 'bad_primes', 'regulator', 'sha']
SEARCH_COLS = ['ainvs', 'lmfdb_label'] + ADDITIONAL_COLS
OUTPUT_COLS = ['lmfdb_label'] + [str(p) for p in primes_first_n(NUM_AP_VALS)] + ADDITIONAL_COLS

MWBSD_COLS = ['lmfdb_label', 'real_period', 'special_value', 'tamagawa_product']

def ap_normalized(E, p):
    ap = E.ap(p)
    normalization_quotient = 2 * p.sqrt()
    return np.float32(round(ap / normalization_quotient, NUM_DECIMAL_PLACES))

# Function to get the data and labels
def get_data(tab):
    data = []
    for curve in tab:
        ainvs = curve['ainvs']
        cond = Integer(curve['conductor'])
        minimal_disc = Integer(curve['absD'])
        bad_primes = curve['bad_primes']
        E = EllipticCurve(ainvs)
        a3 = E.ap(3)  # condition 1b
        condition_1b = (a3 == 0)
        isogeny_degrees = [phi.degree() for phi in E.isogenies_prime_degree()]
        condition_1c = (isogeny_degrees == []) or (isogeny_degrees == [2])

        condition_1d = True  # initialize to True
        for p in bad_primes:
            ord_p_of_min_disc = minimal_disc.valuation(p)
            order_of_order = ord_p_of_min_disc.valuation(p)
            if order_of_order > 0:
                condition_1d = False
                break

        conditions = [condition_1b, condition_1c, condition_1d]

        if all(conditions):
            data.append(curve['lmfdb_label'])

    return data

def foo(cond_bound=20):

    print(f"Generating data file for curves of conductor up to {cond_bound}...")
    output_file = OUTPUT_FILE  # .format(NUM_AP_VALS, 1, MY_LOCAL_LIM-1)
    my_query = {'conductor': {'$lt': cond_bound}, 'semistable' : True}
    the_curves = ecq.search(my_query, projection=SEARCH_COLS, one_per=['lmfdb_iso'])

    data = get_data(the_curves)

    # data = list(the_curves)

    # Export labels to a text file
    with open(output_file, 'w') as f:
        for label in data:
            f.write(f"{label}\n")

    print(f"SUCCESS!!! Data file saved to {output_file}")

print("Working...")

start_time = time.time()

foo(50)

end_time = time.time()
elapsed_time = end_time - start_time

minutes = int(elapsed_time // 60)
seconds = int(elapsed_time % 60)

print(f"Elapsed time: {minutes} minutes {seconds} seconds")
