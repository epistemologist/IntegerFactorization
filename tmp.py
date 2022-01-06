from elliptic_curve import *
from itertools import count

import cProfile

curve = EllipticCurve(2**127 - 1, 12345, 6789)

pts = []

N = 100
for i in count(2**126):
    tmp = curve.get_point(i)
    if tmp:
        x,y = tmp
        pts.append(EllipticCurveAffinePoint(curve, x, y))
        if len(pts) % 1000 == 0: print(len(pts))
    if len(pts) > N: break


