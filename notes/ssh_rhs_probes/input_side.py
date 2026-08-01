"""Is there an INPUT-side property unique to elements 39028/39029/39030?

The exact zeros in u are written by FESOM during step 1; the incoming remapped
restart has none on these three elements.  A fix in the remap tool is only
well-founded if something about the state the remap hands over is distinctive at
exactly these elements.  Test the obvious candidates:

  - did the element exist on the 1905 mesh at all (new element)?
  - did its ulevels / nlevels change at the couple-in?
  - are its three nodes' ulevels / nlevels changed?
  - how shallow is it, and how does that compare to the population?

map_elem.out / map_nod.out carry stable max-mesh IDs across submeshes.
"""
import numpy as np

CO = "/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/couple"
TARGET = [39028, 39029, 39030]


def load(y):
    M = f"{CO}/submesh_{y}-12-31T00:00:00"
    return dict(
        mape=np.loadtxt(f"{M}/map_elem.out").astype(int),
        mapn=np.loadtxt(f"{M}/map_nod.out").astype(int),
        el=np.loadtxt(f"{M}/elem2d.out", skiprows=1).astype(int) - 1,
        ule=np.loadtxt(f"{M}/cavity_elvls.out").astype(int),
        nle=np.loadtxt(f"{M}/elvls.out").astype(int),
        uln=np.loadtxt(f"{M}/cavity_nlvls.out").astype(int),
        nln=np.loadtxt(f"{M}/nlvls.out").astype(int),
        nod=np.loadtxt(f"{M}/nod2d.out", skiprows=1),
    )


a, b = load("1905"), load("1906")
pa = {m: i for i, m in enumerate(a["mape"])}
na = {m: i for i, m in enumerate(a["mapn"])}

print("=== the three elements, 1905 -> 1906 ===")
for e1 in TARGET:
    j = e1 - 1
    mid = b["mape"][j]
    print(f"\n elem {e1} (maxmesh id {mid})  1906: ul={b['ule'][j]} nl={b['nle'][j]}")
    if mid in pa:
        i = pa[mid]
        print(f"   1905: ul={a['ule'][i]} nl={a['nle'][i]}   "
              f"-> dul={b['ule'][j]-a['ule'][i]:+d} dnl={b['nle'][j]-a['nle'][i]:+d}")
    else:
        print("   1905: ELEMENT DID NOT EXIST (new at this couple-in)")
    for k in b["el"][j]:
        m = b["mapn"][k]
        s = f"   node {k+1:7d} (id {m:7d}) 1906 ul={b['uln'][k]} nl={b['nln'][k]}"
        if m in na:
            i = na[m]
            s += f" | 1905 ul={a['uln'][i]} nl={a['nln'][i]} -> dul={b['uln'][k]-a['uln'][i]:+d} dnl={b['nln'][k]-a['nln'][i]:+d}"
        else:
            s += " | NEW NODE"
        print(s)

print("\n=== population context on the 1906 mesh ===")
ne = b["ule"].size
newel = np.array([b["mape"][j] not in pa for j in range(ne)])
dnl = np.array([b["nle"][j] - a["nle"][pa[b["mape"][j]]] if b["mape"][j] in pa else 99
                for j in range(ne)])
dul = np.array([b["ule"][j] - a["ule"][pa[b["mape"][j]]] if b["mape"][j] in pa else 99
                for j in range(ne)])
print(f" new elements at this couple-in     : {int(newel.sum())}")
print(f" elements with nlevels changed      : {int(((dnl != 0) & (dnl != 99)).sum())}")
print(f" elements with ulevels changed      : {int(((dul != 0) & (dul != 99)).sum())}")
wet = b["nle"] - b["ule"]
print(f" wet-level count: min={wet.min()} p1={np.percentile(wet,1):.0f} "
      f"median={np.median(wet):.0f}")
for e1 in TARGET:
    j = e1 - 1
    print(f"   elem {e1}: wet levels={wet[j]}  new={newel[j]}  dnl={dnl[j]}  dul={dul[j]}"
          f"   (shallower than {100*(wet<wet[j]).mean():.2f}% of elements)")

# how many elements are as shallow as 39030, and did any others get u-zeros?
thin = np.where(wet <= 4)[0]
print(f"\n elements with <=4 wet levels: {thin.size}")
print(f"   of these, how many are in the crash trio: {sum(1 for t in TARGET if (t-1) in set(thin.tolist()))}")
