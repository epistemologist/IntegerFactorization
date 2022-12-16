import pstats
t1 = pstats.Stats("times_normal")
t2 = pstats.Stats("times_w_patch")
print(t1.sort_stats("cumtime").print_stats(20))
print(t2.sort_stats("cumtime").print_stats(20))

