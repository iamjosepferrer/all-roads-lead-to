"""Build a compact routing graph of Europe's main roads and ferries for routing in the browser.

Input  (data_raw/): segments.parquet          Overture transportation segments (from extract.py)
                    ne_10m_*.geojson          Natural Earth countries and populated places (download.sh)
                    GHS_POP_E2025_*.tif       GHSL population 2025, 30 arcsec grid (download.sh)
Output (data/):     graph.bin.gz, meta.json   read by index.html

Run from anywhere:  python pipeline/build.py   (about 1 to 2 minutes, ~6 GB RAM peak)
"""
import json, gzip, time, sys, os
from pathlib import Path
from collections import defaultdict
import numpy as np
import pyarrow.parquet as pq
import shapely
from shapely.geometry import shape, Point
from shapely.ops import unary_union
import pyproj
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
RAW = HERE.parent / "data_raw"
OUT = HERE.parent / "data"

T0 = time.time()
def log(*a):
    print(f"[{time.time()-T0:6.1f}s]", *a, flush=True)

SCOPE = set("""AUT BEL BGR HRV CYP CZE DNK EST FIN FRA DEU GRC HUN IRL ITA LVA LTU LUX MLT NLD POL PRT ROU SVK SVN ESP SWE
ISL LIE NOR CHE GBR ALB BIH MNE MKD SRB KOS AND MCO SMR VAT GIB IMN JEY GGY ALD FRO""".split())
BBOX = (-25, 34, 35, 72)
ORIGIN = (4_300_000.0, 3_300_000.0)  # LAEA (EPSG:3035) false origin we subtract, metres

to3035 = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:3035", always_xy=True)

def proj(lon, lat):
    x, y = to3035.transform(lon, lat)
    return (np.asarray(x) - ORIGIN[0]) / 1000.0, (np.asarray(y) - ORIGIN[1]) / 1000.0  # km

# scope polygon
countries = json.load(open(RAW / "ne_10m_admin_0_countries.geojson"))["features"]
scope_polys, context = [], []
bbox_poly = shapely.box(*BBOX)
for f in countries:
    p = f["properties"]
    code = p.get("ADM0_A3")
    g = shape(f["geometry"]).intersection(bbox_poly)
    if g.is_empty:
        continue
    context.append((code, p.get("NAME"), g, code in SCOPE))
    if code in SCOPE:
        scope_polys.append(g)
scope = unary_union(scope_polys)
scope_buf = scope.buffer(0.03)  # ~3 km, keeps border bridges and coastal stubs
shapely.prepare(scope_buf)
log("scope countries", sum(1 for c in context if c[3]))

# segments
tab = pq.read_table(RAW / "segments.parquet")
log("segments", tab.num_rows)
import pyarrow as pa, pyarrow.compute as pc
subtype = np.array(tab.column("subtype").to_pylist(), dtype=object)
klass = np.array([c or "" for c in tab.column("class").to_pylist()], dtype=object)
subcls = np.array([c or "" for c in tab.column("subclass").to_pylist()], dtype=object)
geoms = shapely.from_wkb(tab.column("wkb").to_numpy(zero_copy_only=False))
spd_col = tab.column("spd").combine_chunks()
conn_col = tab.column("connectors").combine_chunks()
c_parent = pc.list_parent_indices(conn_col).to_numpy()
c_flat = pc.list_flatten(conn_col)
c_id = c_flat.field("connector_id").to_numpy(zero_copy_only=False)
c_at = c_flat.field("at").to_numpy(zero_copy_only=False).astype(float)
del tab, conn_col, c_flat

# keep only linestrings
gt = shapely.get_type_id(geoms)
ok = gt == 1
# scope filter: roads by midpoint, ferries need both ends in scope
mid = shapely.line_interpolate_point(geoms, 0.5, normalized=True)
inside_mid = shapely.contains_xy(scope_buf, shapely.get_x(mid), shapely.get_y(mid))
st = shapely.get_point(geoms, 0)
en = shapely.get_point(geoms, -1)
in_st = shapely.contains_xy(scope_buf, shapely.get_x(st), shapely.get_y(st))
in_en = shapely.contains_xy(scope_buf, shapely.get_x(en), shapely.get_y(en))
inside_ends = in_st & in_en
inside_mid = inside_mid | in_st | in_en  # long bridges (Oresund, Great Belt) have their midpoint over water
is_ferry = subtype == "water"
keep = ok & np.where(is_ferry, inside_ends, inside_mid)
idx_keep = np.nonzero(keep)[0]
log("kept segments", len(idx_keep), "ferries", int((is_ferry & keep).sum()))

geoms = geoms[idx_keep]; subtype = None; klass = klass[idx_keep]; subcls = subcls[idx_keep]; is_ferry = is_ferry[idx_keep]
spd = spd_col.take(pa.array(idx_keep)); del spd_col
del mid, st, en, inside_mid, inside_ends, in_st, in_en, gt, ok
newpos = np.full(len(keep), -1, np.int64); newpos[idx_keep] = np.arange(len(idx_keep))
cm = keep[c_parent]
c_seg = newpos[c_parent[cm]]; c_id = c_id[cm]; c_at = np.nan_to_num(c_at[cm])
S = len(geoms)

# vertices (projected) + cumulative fractions
coords, vseg = shapely.get_coordinates(geoms, return_index=True)
vx, vy = proj(coords[:, 0], coords[:, 1])
vx = vx.astype(np.float64); vy = vy.astype(np.float64)
del geoms, coords
same = np.r_[False, vseg[1:] == vseg[:-1]]
step = np.where(same, np.hypot(np.diff(vx, prepend=vx[0]), np.diff(vy, prepend=vy[0])), 0.0)
cum = np.cumsum(step)
seg_start = np.searchsorted(vseg, np.arange(S))
seg_end = np.r_[seg_start[1:], len(vseg)]  # exclusive
cum = cum - cum[seg_start][vseg]
seg_len = cum[seg_end - 1]
frac = np.where(seg_len[vseg] > 0, cum / np.maximum(seg_len[vseg], 1e-12), 0.0)
log("vertices", len(vx))

# connectors -> vertex positions
key = vseg * 2.0 + frac
k = c_seg * 2.0 + np.clip(c_at, 0, 1)
i = np.searchsorted(key, k)
i = np.clip(i, seg_start[c_seg], seg_end[c_seg] - 1)
im = np.maximum(i - 1, seg_start[c_seg])
i = np.where(np.abs(key[im] - k) < np.abs(key[i] - k), im, i)
uid, c_node = np.unique(c_id, return_inverse=True)
cs_seg = c_seg; cs_vert = i; cs_node = c_node
nextid = len(uid)
# ensure nodes at both ends of every segment (vectorised)
v_ends = np.r_[seg_start, seg_end - 1]; s_ends = np.r_[np.arange(S), np.arange(S)]
miss = ~np.isin(v_ends, cs_vert)
mv = v_ends[miss]; ms = s_ends[miss]
rk = np.round(vx[mv], 3) * 1e7 + np.round(vy[mv], 3)
_, inv = np.unique(rk, return_inverse=True)
cs_seg = np.r_[cs_seg, ms].astype(np.int64); cs_vert = np.r_[cs_vert, mv].astype(np.int64)
cs_node = np.r_[cs_node, nextid + inv].astype(np.int64)
N0 = nextid + (inv.max() + 1 if len(inv) else 0)
del v_ends, s_ends, miss, mv, ms, rk, inv, c_id, c_at, c_seg, c_parent
log("nodes (raw)", N0)

# node coordinates
node_x = np.zeros(N0); node_y = np.zeros(N0)
node_x[cs_node] = vx[cs_vert]; node_y[cs_node] = vy[cs_vert]

# speeds (vectorised)
sp_par = pc.list_parent_indices(spd).to_numpy()
sp_flat = pc.list_flatten(spd)
sp_v = sp_flat.field("v").to_numpy(zero_copy_only=False)
sp_u = np.array(sp_flat.field("u").to_pylist(), dtype=object)
sp_b = pc.is_null(sp_flat.field("b")).to_numpy(zero_copy_only=False)
sp_v = np.where(np.isnan(sp_v.astype(float)), 0, sp_v.astype(float))
sp_v = np.where(sp_u == "mph", sp_v * 1.609, sp_v)
sp_v = np.where(sp_b, sp_v, 0)
best = np.zeros(S); np.maximum.at(best, sp_par, sp_v)
del spd, sp_flat, sp_u
base = np.select([klass == "motorway", klass == "trunk", klass == "primary"], [120.0, 100.0, 80.0], 70.0)
base = np.where(subcls == "link", 60.0, base)
seg_speed_arr = np.clip(np.where(best > 0, best, base), 20, 130)
seg_speed_arr = np.where(is_ferry, 35.0, seg_speed_arr)
seg_cls = np.select([is_ferry, klass == "motorway", klass == "trunk"], [0, 1, 2], 3).astype(np.uint8)

# raw edges between consecutive connectors
order = np.lexsort((cs_vert, cs_seg))
cs_seg = cs_seg[order]; cs_vert = cs_vert[order]; cs_node = cs_node[order]
pair = cs_seg[1:] == cs_seg[:-1]
eu = cs_node[:-1][pair]; ev = cs_node[1:][pair]
ea = cs_vert[:-1][pair]; eb = cs_vert[1:][pair]
es = cs_seg[:-1][pair]
elen = cum[eb] - cum[ea]
good = (eu != ev) & (eb > ea)
eu, ev, ea, eb, es, elen = eu[good], ev[good], ea[good], eb[good], es[good], elen[good]
del cs_seg, cs_vert, cs_node, order, pair
straight = np.hypot(node_x[eu] - node_x[ev], node_y[eu] - node_y[ev])
loop = (seg_cls[es] == 0) & (elen > 100) & (elen > 3 * straight)
log("dropped cruise-like ferry edges", int(loop.sum()))
eu, ev, ea, eb, es, elen = eu[~loop], ev[~loop], ea[~loop], eb[~loop], es[~loop], elen[~loop]
etime = elen / seg_speed_arr[es] * 60.0  # minutes
ecls = seg_cls[es]
eferry = ecls == 0
del frac, key, step, vseg, seg_start, seg_end, seg_len
log("raw edges", len(eu))

# ferry boarding penalty: +20 min at every ferry edge end that touches a road
road_nodes = np.zeros(N0, bool); road_nodes[eu[~eferry]] = True; road_nodes[ev[~eferry]] = True
etime = etime + np.where(eferry, 20.0 * (road_nodes[eu].astype(float) + road_nodes[ev].astype(float)), 0.0)
NRAW = len(eu)
extra = dict(u=[], v=[], t=[], l=[], c=[])

def add_extra(u, v, t, l, c):
    for k_, val in zip("utlc", (u, t, l, c)): pass
    extra["u"].append(u); extra["v"].append(v); extra["t"].append(t); extra["l"].append(l); extra["c"].append(c)

def all_edges():
    return (np.r_[eu, np.array(extra["u"], np.int64)], np.r_[ev, np.array(extra["v"], np.int64)])

# connectivity repair
def components(nn):
    U_, V_ = all_edges()
    A = coo_matrix((np.ones(len(U_), np.int8), (U_, V_)), shape=(nn, nn))
    return connected_components(A, directed=False)

used = np.zeros(N0, bool); used[eu] = True; used[ev] = True
ferry_term = np.zeros(N0, bool); ferry_term[eu[eferry]] = True; ferry_term[ev[eferry]] = True
road_only = np.zeros(N0, bool); road_only[eu[~eferry]] = True; road_only[ev[~eferry]] = True
pts_idx = np.nonzero(road_only)[0]
linked = np.zeros(N0, bool)
tree = cKDTree(np.c_[node_x[pts_idx], node_y[pts_idx]])
for it in range(3):
    ncomp, lab = components(N0)
    sizes = np.bincount(lab[used], minlength=ncomp)
    big = int(np.argmax(sizes))
    U_, V_ = all_edges()
    deg = np.bincount(np.r_[U_, V_], minlength=N0)
    dangling_ferry = ferry_term & ~road_only & ~linked
    cand = np.nonzero(used & (((lab != big) & (deg <= 1)) | dangling_ferry))[0]
    R = np.where(ferry_term[cand], 10.0, 8.0)
    d, j = tree.query(np.c_[node_x[cand], node_y[cand]], k=24, distance_upper_bound=10.0, workers=-1)
    best_per_comp = {}
    for row, n in enumerate(cand):
        for dd, jj in zip(d[row], j[row]):
            if not np.isfinite(dd) or dd > R[row]:
                break
            m = pts_idx[jj]
            if lab[m] != lab[n] or dangling_ferry[n]:
                key_ = (lab[n], lab[m]) if not dangling_ferry[n] else ("f", n)
                if key_ not in best_per_comp or dd < best_per_comp[key_][0]:
                    best_per_comp[key_] = (dd, n, m)
                break
    for (dd, n, m) in best_per_comp.values():
        linked[n] = True
        L = dd * 1.3
        add_extra(int(n), int(m), L / 40.0 * 60.0 + 1.0, L, 4)
    log(f"repair pass {it}: components {ncomp}, largest {sizes[big]}, candidates {len(cand)}, links added {len(best_per_comp)}")

# Channel Tunnel car shuttle (Folkestone <-> Coquelles): 35 min crossing + 30 min loading
def nearest_node(lon, lat):
    x, y = proj(lon, lat)
    _, j_ = tree.query([float(x), float(y)])
    return int(pts_idx[j_])
a = nearest_node(1.137, 51.096); b = nearest_node(1.812, 50.935)
add_extra(a, b, 65.0, 50.0, 5)  # class 5 = rail shuttle, kept when ferries are off

# keep largest component
ncomp, lab = components(N0)
sizes = np.bincount(lab[used], minlength=ncomp)
big = int(np.argmax(sizes))
U_, V_ = all_edges()
T_ = np.r_[etime, np.array(extra["t"])]; L_ = np.r_[elen, np.array(extra["l"])]; C_ = np.r_[ecls, np.array(extra["c"], np.uint8)]
keepE = lab[U_] == big
log("largest component nodes", sizes[big], "edges kept", int(keepE.sum()), "of", len(keepE))
rid = np.nonzero(keepE)[0]  # raw edge id (>= NRAW means extra straight link)

# drop parallel duplicates (keep fastest)
pk = np.minimum(U_[rid], V_[rid]).astype(np.int64) * N0 + np.maximum(U_[rid], V_[rid])
o = np.lexsort((T_[rid], pk))
first = np.r_[True, pk[o][1:] != pk[o][:-1]]
rid = rid[o[first]]
eu2, ev2, et, el, ec = U_[rid], V_[rid], T_[rid], L_[rid], C_[rid]
log("edges after dedup", len(rid))

# contract degree-2 chains (numba)
from numba import njit

@njit(cache=True)
def contract(N0, eu2, ev2, et, el, ec, adj_off, adj_e, is_junction):
    M = len(eu2)
    visited = np.zeros(M, np.bool_)
    out_u = np.empty(M, np.int64); out_v = np.empty(M, np.int64)
    out_t = np.empty(M, np.float64); out_l = np.empty(M, np.float64); out_c = np.empty(M, np.uint8)
    part_e = np.empty(M, np.int64); part_f = np.empty(M, np.bool_); part_off = np.zeros(M + 1, np.int64)
    ne = 0; np_ = 0
    for n in range(N0):
        if not is_junction[n] or adj_off[n] == adj_off[n + 1]:
            continue
        for jj in range(adj_off[n], adj_off[n + 1]):
            i = adj_e[jj]
            if visited[i]:
                continue
            cur = n; t = 0.0; l = 0.0; cl = 9; anyroad = False
            while True:
                visited[i] = True
                fwd = eu2[i] == cur
                nxt = ev2[i] if fwd else eu2[i]
                part_e[np_] = i; part_f[np_] = fwd; np_ += 1
                t += et[i]; l += el[i]
                if ec[i] != 0:
                    anyroad = True
                    if ec[i] < cl:
                        cl = ec[i]
                cur = nxt
                if is_junction[cur] or cur == n:
                    break
                a0 = adj_e[adj_off[cur]]; a1 = adj_e[adj_off[cur] + 1]
                i = a1 if a0 == i else a0
                if visited[i]:
                    break
            out_u[ne] = n; out_v[ne] = cur; out_t[ne] = t; out_l[ne] = l
            out_c[ne] = cl if anyroad else 0
            ne += 1; part_off[ne] = np_
    return out_u[:ne], out_v[:ne], out_t[:ne], out_l[:ne], out_c[:ne], part_e[:np_], part_f[:np_], part_off[:ne + 1]

@njit(cache=True)
def assemble(part_e, part_f, part_off, raw_a, raw_b, raw_u, raw_v, vx, vy, node_x, node_y):
    ne = len(part_off) - 1
    # count vertices
    total = 0
    for k in range(ne):
        for p in range(part_off[k], part_off[k + 1]):
            r = part_e[p]
            nv = raw_b[r] - raw_a[r] + 1 if raw_a[r] >= 0 else 2
            total += nv if p == part_off[k] else nv - 1
    xy = np.empty((total, 2), np.float64); idx = np.empty(total, np.int64)
    w = 0
    for k in range(ne):
        for p in range(part_off[k], part_off[k + 1]):
            r = part_e[p]; fwd = part_f[p]; skip = p != part_off[k]
            if raw_a[r] >= 0:
                a = raw_a[r]; b = raw_b[r]; nv = b - a + 1
                for q in range(nv):
                    if skip and q == 0:
                        continue
                    s = a + q if fwd else b - q
                    xy[w, 0] = vx[s]; xy[w, 1] = vy[s]; idx[w] = k; w += 1
            else:
                u = raw_u[r]; v = raw_v[r]
                if not fwd:
                    u, v = v, u
                if not skip:
                    xy[w, 0] = node_x[u]; xy[w, 1] = node_y[u]; idx[w] = k; w += 1
                xy[w, 0] = node_x[v]; xy[w, 1] = node_y[v]; idx[w] = k; w += 1
    return xy[:w], idx[:w]

for it in range(6):
    deg = np.bincount(np.r_[eu2, ev2], minlength=N0)
    leaf = (deg[eu2] == 1) | (deg[ev2] == 1)
    drop = leaf & (el < 2.0) & (ec != 0)
    if not drop.any():
        break
    k_ = ~drop
    rid, eu2, ev2, et, el, ec = rid[k_], eu2[k_], ev2[k_], et[k_], el[k_], ec[k_]
    log(f"prune pass {it}: dropped {int(drop.sum())} stubs")
M = len(rid)
deg = np.bincount(np.r_[eu2, ev2], minlength=N0)
ends = np.r_[eu2, ev2]; eids = np.r_[np.arange(M), np.arange(M)]
o = np.argsort(ends, kind="stable")
adj_e = eids[o]; adj_off = np.r_[0, np.cumsum(np.bincount(ends, minlength=N0))]
del ends, eids, o
fer_ = ec == 0
touch_f = np.zeros(N0, bool); touch_f[eu2[fer_]] = True; touch_f[ev2[fer_]] = True
touch_r = np.zeros(N0, bool); touch_r[eu2[~fer_]] = True; touch_r[ev2[~fer_]] = True
junc = (deg != 2) | (touch_f & touch_r)  # never merge a ferry with a road link
cu, cv, ct, cl_, cc, part_e, part_f, part_off = contract(N0, eu2.astype(np.int64), ev2.astype(np.int64), et.astype(np.float64),
                                                       el.astype(np.float64), ec.astype(np.uint8), adj_off, adj_e, junc)
newE = list(zip(cu.tolist(), cv.tolist(), ct.tolist(), cl_.tolist(), cc.tolist()))
log("contracted edges", len(newE))

used_nodes = np.unique(np.r_[cu, cv])
remap = {int(n): i for i, n in enumerate(used_nodes)}
NN = len(used_nodes)
nx = node_x[used_nodes].astype(np.float32); ny = node_y[used_nodes].astype(np.float32)

nextra = len(extra["u"])
raw_a = np.r_[ea, np.full(nextra, -1)][rid].astype(np.int64)
raw_b = np.r_[eb, np.full(nextra, -1)][rid].astype(np.int64)
xy, gidx = assemble(part_e, part_f, part_off, raw_a, raw_b, U_[rid].astype(np.int64), V_[rid].astype(np.int64), vx, vy, node_x, node_y)
del vx, vy, cum
log("assembled vertices", len(xy))
lines = shapely.linestrings(xy, indices=gidx)
del xy, gidx
simp = shapely.simplify(lines, 0.2, preserve_topology=False)
del lines
def _spread(v):
    v = v.astype(np.uint64) & np.uint64(0xFFFFF)
    out = np.zeros_like(v)
    for b in range(20):
        out |= ((v >> np.uint64(b)) & np.uint64(1)) << np.uint64(2 * b)
    return out
qnx = np.round(nx.astype(np.float64) * 20).astype(np.int64); qny = np.round(ny.astype(np.float64) * 20).astype(np.int64)
mort = _spread(qnx - qnx.min()) | (_spread(qny - qny.min()) << np.uint64(1))
nord = np.argsort(mort, kind="stable")
rank = np.empty(NN, np.int64); rank[nord] = np.arange(NN)
nx, ny, qnx, qny = nx[nord], ny[nord], qnx[nord], qny[nord]
Ui = rank[np.searchsorted(used_nodes, cu)]; Vi = rank[np.searchsorted(used_nodes, cv)]
perm = np.lexsort((Vi, Ui))
Ui, Vi = Ui[perm], Vi[perm]
ct, cl_, cc = ct[perm], cl_[perm], cc[perm]
simp = simp[perm]
newE = list(zip(Ui.tolist(), Vi.tolist(), ct.tolist(), cl_.tolist(), cc.tolist()))
remap = {i: i for i in range(NN)}
pcoords, pidx = shapely.get_coordinates(simp, return_index=True)
counts = np.bincount(pidx, minlength=len(newE))
log("nodes", NN, "edges", len(newE), "path vertices", len(pcoords))

# population per edge (GHSL 2025)
import rasterio
from rasterio.windows import from_bounds
from rasterio import features
with rasterio.open(RAW / "GHS_POP_E2025_GLOBE_R2023A_4326_30ss_V1_0.tif") as src:
    win = from_bounds(*BBOX, transform=src.transform)
    arr = src.read(1, window=win).astype(np.float64)
    tr = src.window_transform(win)
arr[~np.isfinite(arr) | (arr < 0)] = 0
mask = features.geometry_mask([scope], out_shape=arr.shape, transform=tr, invert=True)
arr[~mask] = 0
prow, pcol = np.nonzero(arr > 0)
plon = tr.c + (pcol + 0.5) * tr.a
plat = tr.f + (prow + 0.5) * tr.e
pval = arr[prow, pcol]
px, py = proj(plon, plat)
# assign each pixel to nearest path vertex -> its edge
vtree = cKDTree(pcoords)
_, vi = vtree.query(np.c_[px, py], workers=-1)
epop = np.bincount(pidx[vi], weights=pval, minlength=len(newE))
log("population assigned", f"{pval.sum()/1e6:.1f} M people")

# cities
places = json.load(open(RAW / "ne_10m_populated_places_simple.geojson"))["features"]
ntree = cKDTree(np.c_[nx, ny])
cities = []
for f in places:
    p = f["properties"]; lon, lat = f["geometry"]["coordinates"]
    if not (BBOX[0] < lon < BBOX[2] and BBOX[1] < lat < BBOX[3]):
        continue
    if not scope.contains(Point(lon, lat)):
        continue
    pop = p.get("pop_max") or 0
    cap = p.get("adm0cap") == 1 or p.get("featurecla", "").startswith("Admin-0 capital")
    if pop < 40000 and not cap:
        continue
    x, y = proj(lon, lat)
    d, j = ntree.query([float(x), float(y)])
    if d > 25:
        continue
    cities.append(dict(name=p["name"], cc=p.get("adm0_a3") or p.get("sov_a3"), pop=int(pop),
                       cap=bool(cap), n=int(j), x=round(float(x), 2), y=round(float(y), 2)))
cities.sort(key=lambda c: -c["pop"])
log("cities", len(cities))

# land outlines for the basemap
land = []
for code, name, g, ins in context:
    g2 = shapely.simplify(shapely.transform(g, lambda c: np.c_[proj(c[:, 0], c[:, 1])]), 1.5)
    polys = [g2] if g2.geom_type == "Polygon" else list(getattr(g2, "geoms", []))
    rings = []
    for pg in polys:
        if pg.geom_type != "Polygon" or pg.area < 30:
            continue
        rings.append([[round(a, 1), round(b, 1)] for a, b in pg.exterior.coords])
    if rings:
        land.append(dict(cc=code, name=name, in_scope=ins, rings=rings))

# write binary (byte-shuffled columns, gzip)
import os
os.makedirs(OUT, exist_ok=True)
E2 = len(newE)
def shuffle(arr):
    a8 = np.ascontiguousarray(arr).view(np.uint8).reshape(-1, arr.dtype.itemsize)
    return a8.T.copy().tobytes()
q = np.round(pcoords * 20).astype(np.int64)  # 50 m units
counts = counts.astype(np.int64)
starts = np.r_[0, np.cumsum(counts)[:-1]]
d = np.diff(q, axis=0, prepend=q[:1])
d[starts, 0] = q[starts, 0] - qnx[Ui]; d[starts, 1] = q[starts, 1] - qny[Ui]
assert np.abs(d).max() < 32767, np.abs(d).max()
assert counts.max() < 65535
dU = np.diff(Ui, prepend=0).astype(np.uint32)
dV = (Vi - Ui).astype(np.int32)
TT = np.clip(np.round(ct * 10), 0, 65535).astype(np.uint16)      # 0.1 min
LL = np.clip(np.round(cl_ * 10), 0, 65535).astype(np.uint16)     # 0.1 km
PP = np.round(epop).astype(np.uint32)
dqx = np.diff(qnx, prepend=0).astype(np.int32); dqy = np.diff(qny, prepend=0).astype(np.int32)
header = np.array([1, NN, E2, len(pcoords)], np.uint32)
parts = [header.tobytes(), shuffle(dqx), shuffle(dqy), shuffle(dU), shuffle(dV), shuffle(TT), shuffle(LL), shuffle(PP),
         cc.astype(np.uint8).tobytes(), shuffle(counts.astype(np.uint16)), shuffle(d[:, 0].astype(np.int16).copy()),
         shuffle(d[:, 1].astype(np.int16).copy())]
blob = b"".join(parts)
for f_ in os.listdir(OUT):
    if f_.startswith("graph"):
        os.remove(OUT / f_)
with gzip.open(OUT / "graph.bin.gz", "wb", compresslevel=9) as fh:
    fh.write(blob)
log("graph raw MB", len(blob) / 1e6, "gz MB", os.path.getsize(OUT / "graph.bin.gz") / 1e6)
meta = dict(origin=ORIGIN, crs="EPSG:3035", units="km", quant=0.05,
            built=time.strftime("%Y-%m-%d"), source="Overture Maps 2026-10-21.0 (OpenStreetMap, ODbL); GHSL POP 2025",
            cities=cities, land=land)
json.dump(meta, open(OUT / "meta.json", "w"), separators=(",", ":"))
log("meta MB", os.path.getsize(OUT / "meta.json") / 1e6)
log("done")
