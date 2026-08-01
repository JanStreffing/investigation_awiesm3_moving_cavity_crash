"""The root cause, as a pure mesh predicate.

fer_solve_Gamma (oce_fer_gm.F90) solves the GM streamfunction over
    nzmin = ulevels_nod2D_max(n) = max over surrounding elements of ulevels
    nzmax = nlevels_nod2D_min(n) = min over surrounding elements of nlevels
and then unconditionally executes

    tr(:,nzmax) = tp(:,nzmax)

If nzmin > nzmax the range is inverted: both DO loops run zero times, tp(:,nzmax)
is never assigned, and an uninitialised local automatic array element is copied
straight into fer_gamma(:,nzmax,n). fer_gamma2vel then divides gamma differences
by helem to form fer_uv, solve_tracers_ale adds fer_uv to UV, advects tracers with
it, and subtracts it back -- leaving temp ~1e250 and, via (x+y)-y, exact zeros in u.

Predicate: max(ulevels(elems at n)) > min(nlevels(elems at n)).
It needs a node that simultaneously touches a deep-cavity element and a very
shallow one - a cavity front sitting on steep bathymetry.
"""
import numpy as np

CO = "/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/couple"


def check(tag, mdir):
    elem = np.loadtxt(f"{mdir}/elem2d.out", skiprows=1).astype(int) - 1
    ule = np.loadtxt(f"{mdir}/cavity_elvls.out").astype(int)   # ulevels(elem)
    nle = np.loadtxt(f"{mdir}/elvls.out").astype(int)          # nlevels(elem)
    nod = np.loadtxt(f"{mdir}/nod2d.out", skiprows=1)
    N = nod.shape[0]

    umax = np.zeros(N, int)
    nmin = np.full(N, 10**6, int)
    for e in range(elem.shape[0]):
        u, b = ule[e], nle[e]
        for k in elem[e]:
            if u > umax[k]:
                umax[k] = u
            if b < nmin[k]:
                nmin[k] = b
    touched = nmin < 10**6
    bad = np.where(touched & (umax > nmin))[0]
    print(f"\n=== {tag} ===  nodes={N}")
    print(f"  nodes with ulevels_nod2D_max > nlevels_nod2D_min : {bad.size}")
    for n in bad:
        et = np.where(np.any(elem == n, axis=1))[0]
        print(f"    node {n+1:7d} lon {nod[n,1]:8.3f} lat {nod[n,2]:8.3f}  "
              f"nzmin={umax[n]} nzmax={nmin[n]}  (INVERTED by {umax[n]-nmin[n]})")
        print(f"      elements {list(et+1)}  ulevels {list(ule[et])}  nlevels {list(nle[et])}")
    # near misses: equal is degenerate too (single-level solve, but tp(:,nzmax) IS set)
    eq = int(np.sum(touched & (umax == nmin)))
    print(f"  nodes with ulevels_nod2D_max == nlevels_nod2D_min : {eq} (degenerate but safe)")
    return bad.size


tot = 0
for tag in ("1901", "1902", "1903", "1904", "1905", "1906"):
    tot += check(f"submesh {tag}", f"{CO}/submesh_{tag}-12-31T00:00:00")
check("base mesh (pre-coupling)", "/work/ab0246/a270092/input//fesom2/pism_cavity_ini/submesh/")
print(f"\ntotal inverted-range nodes across the six coupled submeshes: {tot}")
