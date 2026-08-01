import numpy as np, os, glob
os.environ.setdefault("HDF5_USE_FILE_LOCKING","FALSE")
import netCDF4 as nc

MESH="/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/couple/submesh_1906-12-31T00:00:00"
R="/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/run_awiesm3_19060101-19061231/work/fesom.1905.oce.restart"
LON0,LAT0=-100.382545,-72.8639145

# --- mesh ---
nod=np.loadtxt(os.path.join(MESH,"nod2d.out"),skiprows=1)
# cols: idx, lon, lat, flag
mlon=nod[:,1]; mlat=nod[:,2]; N=nod.shape[0]
elem=np.loadtxt(os.path.join(MESH,"elem2d.out"),skiprows=1).astype(int)-1  # 0-based node idx
print("mesh nodes",N,"elems",elem.shape[0])

# nearest node to crash lon/lat
d=(mlon-LON0)**2+(mlat-LAT0)**2
i0=int(np.argmin(d))
print(f"crash node global idx={i0+1} lon={mlon[i0]:.4f} lat={mlat[i0]:.4f} dist={d[i0]**.5:.4f}deg")

# neighbor nodes within 30 km (~0.3 deg lon at this lat, use 0.5 deg box)
sel=np.where((np.abs(mlon-LON0)<0.6)&(np.abs(mlat-LAT0)<0.25))[0]
print("nodes in box:",sel.size)

def rd(f):
    p=os.path.join(R,f)
    if not os.path.exists(p): return None
    with nc.Dataset(p) as ds:
        v=list(ds.variables)
        # pick the data var (last, non-dim)
        name=[x for x in v if x not in ds.dimensions][-1]
        arr=ds.variables[name][:]
        return name,np.asarray(arr)

for f in ["ssh.nc","hbar.nc","u.nc","v.nc","urhs_AB.nc","vrhs_AB.nc","w.nc"]:
    r=rd(f)
    if r is None:
        print(f,"MISSING"); continue
    nm,a=r
    print(f"\n=== {f} var={nm} shape={a.shape} ===")
    a=np.squeeze(a)
    # node field?
    if a.ndim==1 and a.shape[0]==N:
        print("  NODE field. crash node value=",a[i0])
        print("  box nodes min/max/mean:",np.nanmin(a[sel]),np.nanmax(a[sel]),np.nanmean(a[sel]))
        print("  GLOBAL min/max:",np.nanmin(a),np.nanmax(a))
    elif a.ndim==2 and a.shape[1]==N:
        col=a[:,i0]
        print("  NODE 3D field. crash node col[:6]=",col[:6])
        print("  box surface min/max:",np.nanmin(a[0,sel]),np.nanmax(a[0,sel]))
        print("  GLOBAL abs max:",np.nanmax(np.abs(a)))
    else:
        # element field
        ne=elem.shape[0]
        if a.ndim==2 and a.shape[1]==ne:
            # elements touching crash node / box nodes
            emask=np.any(np.isin(elem,sel),axis=1)
            print("  ELEM 3D field. elems touching box:",emask.sum())
            sub=a[:,emask]
            print("  box elems abs max:",np.nanmax(np.abs(sub)))
            # per-level max in box
            lvlmax=np.nanmax(np.abs(sub),axis=1)
            print("  box elems |val| per-level[:18]:",np.round(lvlmax[:18],3))
            print("  GLOBAL abs max:",np.nanmax(np.abs(a)), "at elem",np.unravel_index(np.nanargmax(np.abs(a)),a.shape))
        elif a.ndim==1 and a.shape[0]==ne:
            emask=np.any(np.isin(elem,sel),axis=1)
            print("  ELEM 2D field. box elems abs max:",np.nanmax(np.abs(a[emask])),"global",np.nanmax(np.abs(a)))
        else:
            print("  UNRECOGNIZED shape vs N=",N,"ne=",ne)
