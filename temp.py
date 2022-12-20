from multiprocessing import Pool
from itertools import count
from random import random


def gen_args():
    for i in count(1):
        yield i

def f(x):
    if random() < 0.3:
        return [x**2]
    else:
        return None

out = []
with Pool(processes=2) as p:
    for result in p.imap_unordered(f, gen_args()):
        print(len(out))
        if result:
            out.extend(result)
        if len(out) > 1000:
            p.terminate()
            break

print(out)

