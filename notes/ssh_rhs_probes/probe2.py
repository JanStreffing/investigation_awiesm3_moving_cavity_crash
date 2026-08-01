import numpy as np, os, netCDF4 as nc
os.environ["HDF5_USE_FILE_LOCKING"]="FALSE"
MESH="/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/couple/submesh_1906-12-31T00:00:00"
R="/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/run_awiesm3_19060101-19061231/work/fesom.1905.oce.restart"
LON0,LAT0=-100.382545,-72.8639145
nod=np.loadtxt(os.path.join(MESH,"nod2d.out"),skiprows=1)
mlon,mlat=nod[:,1],nod[:,2]; N=nod.shape[0]
elem=np.loadtxt(os.path.join(MESH,"elem2d.out"),skiprows=1).astype(int)-1
# element area on sphere (approx, planar in deg*cos scaled to m)
R_E=6371000.0; rad=np.pi/180
x=R_E*rad*mlon*np.cos(LAT0*rad); y=R_E*rad*mlat
def area(e):
    x1,y1=x[e[:,0]],y[e[:,0]]; x2,y2=x[e[:,1]],y[e[:,1]]; x3,y3=x[e[:,2]],y[e[:,2]]
    return 0.5*np.abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))
ea=area(elem)
i0=int(np.argmin((mlon-LON0)**2+(mlat-LAT0)**2))
# elements touching crash node
etouch=np.where(np.any(elem==i0,axis=1))[0]
# box
sel=np.where((np.abs(mlon-LON0)<0.6)&(np.abs(mlat-LAT0)<0.25))[0]
ebox=np.where(np.any(np.isin(elem,sel),axis=1))[0]
def rd(f,var):
    with nc.Dataset(os.path.join(R,f)) as ds: return np.squeeze(np.asarray(ds.variables[var][:]))
ur=rd("urhs_AB.nc","urhs_AB"); vr=rd("vrhs_AB.nc","vrhs_AB")
ssh=rd("ssh.nc","ssh"); uu=rd("u.nc","u")
print("crash node idx(1-based)=",i0+1,"lon %.3f lat %.3f"%(mlon[i0],mlat[i0]))
print("\n--- box NODE ssh (see fresh-node jumps) ---")
order=sel[np.argsort(ssh[sel])]
for n in order:
    print(f"  node {n+1:7d} lon {mlon[n]:8.3f} lat {mlat[n]:8.3f}  ssh={ssh[n]:8.4f}  flag={int(nod[n,3])}")
print("\n--- crash-node elements: area, surface urhs/vrhs, ratio to neighbors ---")
print("elem_area range in box: %.3e .. %.3e  (ratio %.1f)"%(ea[ebox].min(),ea[ebox].max(),ea[ebox].max()/ea[ebox].min()))
for e in etouch:
    nds=elem[e]
    print(f"  elem {e+1:7d} area={ea[e]:.3e}  urhs0={ur[0,e]:10.2f} vrhs0={vr[0,e]:10.2f}  |u0|={abs(uu[0,e]):.4f}")
print("\n--- surface urhs_AB / area  (accel-like) across box elems: is it spatially jagged? ---")
acc=ur[0,ebox]/ea[ebox]
print("  urhs0/area  min/max/mean:",acc.min(),acc.max(),acc.mean())
print("  raw urhs0   min/max     :",ur[0,ebox].min(),ur[0,ebox].max())
# global: fraction of elems with |urhs0/area| comparable
gacc=np.abs(ur[0,:]/ea)
print("  global |urhs0/area| 99.9pct / max:",np.nanpercentile(gacc,99.9),np.nanmax(gacc))
print("  crash elems |urhs0/area|:",np.round(np.abs(ur[0,etouch]/ea[etouch]),4))
