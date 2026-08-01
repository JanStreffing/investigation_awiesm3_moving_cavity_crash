"""Are the regenerated halo-communication arrays self-consistent?

exchange_elem(UV) at the end of update_vel is the only thing that writes UV after
the momentum arithmetic has produced it, and the crash cluster straddles ranks
1603/1606 - the two that reported the blow-up. dist_NNNN, including these com_info
arrays, is rebuilt by the moving-cavity workflow every cycle, which is consistent
with a failure mode absent from a decade of fixed-mesh FESOM.

Invariant: for every ordered pair (A,B), the number of items A sends to B must equal
the number B receives from A, for both the node and the element exchange. A mismatch
means a rank writes or reads the wrong slice of the halo buffer.

com_info<rank>.out layout, tokens in order:
   rank
   rPEnum, rPE(rPEnum), rptr(rPEnum+1), rlist(rptr(end)-1)
   sPEnum, sPE(sPEnum), sptr(sPEnum+1), slist(sptr(end)-1)
   ... repeated for elements
"""
import glob, numpy as np

CO = "/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/couple"


class Toks:
    def __init__(self, path):
        self.v = np.fromstring(open(path).read(), dtype=np.int64, sep=' ')
        self.i = 0

    def take(self, n=1):
        out = self.v[self.i:self.i + n]
        self.i += n
        return out

    def left(self):
        return self.v.size - self.i


def block(t):
    """one recv-or-send block -> {peer: count}"""
    n = int(t.take()[0])
    pes = t.take(n)
    ptr = t.take(n + 1)
    t.take(int(ptr[-1]) - 1)          # the index list itself
    return {int(pes[k]): int(ptr[k + 1] - ptr[k]) for k in range(n)}


def parse(path):
    t = Toks(path)
    rank = int(t.take()[0])
    out = {}
    for key in ("node_recv", "node_send", "elem_recv", "elem_send"):
        if t.left() <= 0:
            break
        out[key] = block(t)
    return rank, out


for tag in ("1901", "1902", "1903", "1904", "1905", "1906"):
    dists = glob.glob(f"{CO}/submesh_{tag}-12-31T00:00:00/dist_*")
    if not dists:
        continue
    files = sorted(glob.glob(f"{dists[0]}/com_info*.out"))
    if not files:
        continue
    info = {}
    for f in files:
        try:
            r, o = parse(f)
            info[r] = o
        except Exception as e:
            print(f"  parse failure {f}: {e}")
    bad_n, bad_e, missing = 0, 0, 0
    for a, oa in info.items():
        for kind, (skey, rkey) in (("node", ("node_send", "node_recv")),
                                   ("elem", ("elem_send", "elem_recv"))):
            for b, cnt in oa.get(skey, {}).items():
                ob = info.get(b)
                if ob is None or rkey not in ob:
                    missing += 1
                    continue
                if ob[rkey].get(a, -1) != cnt:
                    if kind == "node":
                        bad_n += 1
                    else:
                        bad_e += 1
    print(f"submesh {tag}: {len(info)} ranks   send/recv count mismatches: "
          f"nodes={bad_n}  elements={bad_e}   unresolved peers={missing}")
