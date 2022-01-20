"""
from collections import defaultdict
from utility import get_primes
from tqdm import tqdm

B1 = 10000
B2 = 50000
D = 100

primes = get_primes(B2 + 2*D)
prime_map = {r: [q for q in primes if r+2 <= q <= r+2*D] for r in range(B1-1, B2, 2*D) }

r_s = iter(range(B1-1, B2, 2*D))
prime_map2 = defaultdict(list)

r = next(r_s)

for p in tqdm(primes):
    if r+2 <= p <= r+2*D:
        prime_map2[r].append(p)
    else:
        r = next(r_s)
"""

from utility import tonelli_shanks

def sqrt_mod_p_squared(n, p):
    # Based on https://mathoverflow.net/a/223806
    x = tonelli_shanks(n, p)
    if x is None:
        return None
    roots = []
    for root in [x, -x%p]:
        roots.append((pow(root, p, p*p) * pow(n, (p*p-p+1) // 2, p*p)) % (p*p))
    return roots

