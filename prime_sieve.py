from typing import List
from tqdm import tqdm
from timeit import timeit
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from multiprocessing import cpu_count
from itertools import chain
try:
    from math import isqrt
except:
    isqrt = lambda x: int(x**0.5)

bit_mask = int

# Utility functions

def naive_sieve_of_eratosthenes(n: int) -> List[int]:
    is_prime = [True] * (n-2)
    for i in range(2, int(n**0.5)):
        if is_prime[i-2]:
            for j in range(i*i, n, i):
                is_prime[j-2] = False
    primes = []
    for i in range(len(is_prime)):
        if is_prime[i]: primes.append(i+2)
    return primes

def gen_bit_mask(length: int, spacing: int) -> bit_mask:
    mask = 0
    for i in range(spacing-1, length+spacing, spacing):
        mask |= (1 << i)
    return bit_mask(mask)

def get_small_primes(n: int) -> List[int]:
    small_primes = naive_sieve_of_eratosthenes(100)
    if n < 100:
        return [p for p in small_primes if p <= n]
    else:
        return naive_sieve_of_eratosthenes(n)

class FastPrimeSieve:
    def __init__(self, n: int, threading: bool = False, progress_bar: bool = True):
        self.n = n
        self.threading = threading
        self.chunk_size = isqrt(n)
        self.small_primes = get_small_primes(self.chunk_size)
        self.bit_masks = {p: gen_bit_mask(self.chunk_size, p) for p in self.small_primes}
        self.primes = list(self.small_primes)
        self.full_chunk = bit_mask(2**self.chunk_size - 1)
        # Progress bar is optional to avoid dependencies (this can run on pure Python)
        self.progress_bar = tqdm(total=len(range(1 + self.chunk_size, self.n, self.chunk_size))) if progress_bar else None
    def get_bit_mask(self, chunk_start: int, p: int) -> int:
        inv_bit_mask = self.bit_masks[p] >> (p - 1 - chunk_start)
        return inv_bit_mask ^ self.full_chunk
    def get_primes(self) -> None:
        if self.threading:
            with ThreadPoolExecutor(max_workers=cpu_count()) as executor:
                self.primes = list(chain.from_iterable(executor.map(self.process_chunk_threaded, range(1 + self.chunk_size, self.n, self.chunk_size))))
            print("Sorting primes...")
            self.primes.sort()
            print("Done!")
        else:
            for chunk_start in range(1 + self.chunk_size, self.n, self.chunk_size):
                self.primes.extend(self.process_chunk(chunk_start))
                if self.progress_bar: self.progress_bar.update()
    def process_chunk(self, chunk_start: int) -> List[int]:
        curr_chunk = self.full_chunk
        primes = []
        for p in self.small_primes:
            curr_chunk &= self.get_bit_mask(-chunk_start % p, p)
        for i in range(self.chunk_size):
            if curr_chunk & (1 << i):
                primes.append(chunk_start + i)
        return primes
    def process_chunk_threaded(self, chunk_start: int) -> List[int]:
        result = self.process_chunk(chunk_start)
        if self.progress_bar: self.progress_bar.update()
        return result

