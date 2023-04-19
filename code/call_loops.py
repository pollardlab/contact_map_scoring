#!/usr/env/bin python

import gc
import numpy as np
from collections import defaultdict
from peakachu import peakacluster


def call_loops(score, ther=0.6):
    X = score[score[:, 2] >= ther, :]
    r = X[:, 0].astype(int)
    c = X[:, 1].astype(int)
    p = X[:, 2].astype(float)
    raw = X[:, 3].astype(float)
    d = c - r
    tmpr, tmpc, tmpp, tmpraw, tmpd = r, c, p, raw, d
    matrix = {(r[i], c[i]): p[i] for i in range(len(r))}
    count = 40001
    while count > 40000:
        D = defaultdict(float)
        P = defaultdict(float)
        unique_d = list(set(tmpd.tolist()))
        for distance in unique_d:
            dx = (tmpd == distance)
            dr, dc, dp, draw = tmpr[dx], tmpc[dx], tmpp[dx], tmpraw[dx]
            dx = (dp > np.percentile(dp, 10))
            dr, dc, dp, draw = dr[dx], dc[dx], dp[dx], draw[dx]
            for i in range(dr.size):
                D[(dr[i], dc[i])] += draw[i]
                P[(dr[i], dc[i])] += dp[i]
        count = len(D.keys())
        tmpr = np.array([i[0] for i in P.keys()])
        tmpc = np.array([i[1] for i in P.keys()])
        tmpp = np.array([P.get(i) for i in P.keys()])
        tmpraw = np.array([D.get(i) for i in P.keys()])
        tmpd = tmpc - tmpr

    del X
    gc.collect()
    final_list = peakacluster.local_clustering(D, 1)
    final_list = [i[0] for i in final_list]
    r = [i[0] for i in final_list]
    c = [i[1] for i in final_list]
    p = np.array([matrix.get((r[i], c[i])) for i in range(len(r))])
    if len(r) > 7000:
        sorted_index = np.argsort(p)
        r = [r[i] for i in sorted_index[-7000:]]
        c = [c[i] for i in sorted_index[-7000:]]
    loops = []
    for i in range(len(r)):
        P = matrix.get((r[i], c[i]))
        loop = str(r[i]) + "," + str(c[i]) + ":" + str(P)
        loops.append(loop)
    loops_list = ";".join(loops)
    return r, c, loops_list
