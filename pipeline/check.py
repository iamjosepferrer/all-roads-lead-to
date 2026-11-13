"""Sanity check the built graph: drive times from Brussels, with and without ferries.

Run after build.py:  python pipeline/check.py
"""
import gzip, json
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra

DATA = Path(__file__).resolve().parent.parent / "data"
b = gzip.open(DATA / "graph.bin.gz").read()
_, NN, E, PV = np.frombuffer(b[:16], np.uint32)
off = 16

def take(n, dt):
    """Read one byte-shuffled column."""
    global off
    s = np.dtype(dt).itemsize
    a = np.frombuffer(b, np.uint8, n * s, off).reshape(s, n).T.copy().view(dt).ravel()
    off += n * s
    return a

take(NN, np.int32); take(NN, np.int32)  # node coordinates, not needed here
U = np.cumsum(take(E, np.uint32).astype(np.int64)); V = U + take(E, np.int32)
T = take(E, np.uint16) / 10; take(E, np.uint16); P = take(E, np.uint32)
C = np.frombuffer(b, np.uint8, E, off)

meta = json.load(open(DATA / "meta.json"))
cities = {c["name"]: c["n"] for c in meta["cities"]}
src = cities["Brussels"]

def times(mask):
    A = coo_matrix((T[mask] + 1e-6, (U[mask], V[mask])), shape=(NN, NN)).tocsr()
    return dijkstra(A, directed=False, indices=src)

d_all = times(np.ones(E, bool))
d_road = times(C != 0)  # class 0 = ferry
print(f"{NN:,} junctions, {E:,} links, {P.sum() / 1e6:.0f} M people")
print(f"{'from Brussels':16s} {'ferries':>8s} {'roads only':>11s}")
for name in ["London", "Amsterdam", "Warsaw", "Lisbon", "Stockholm", "Oslo", "Helsinki", "Athens", "Dublin", "Palermo"]:
    if name in cities:
        n = cities[name]
        road = f"{d_road[n] / 60:.1f} h" if np.isfinite(d_road[n]) else "no road"
        print(f"{name:16s} {d_all[n] / 60:7.1f} h {road:>11s}")
reach = np.isfinite(d_road[U]) & np.isfinite(d_road[V])
print(f"people cut off without ferries: {100 * P[~reach].sum() / P.sum():.1f}%")
