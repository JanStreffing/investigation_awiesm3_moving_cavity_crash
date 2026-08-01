import numpy as np, os, netCDF4 as nc
os.environ["HDF5_USE_FILE_LOCKING"]="FALSE"
CO="/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/couple"
R="/work/ab0246/a270092/runtime/awiesm3-develop-is/fgdbg02/run_awiesm3_19060101-19061231/work/fesom.1905.oce.restart"
M=f"{CO}/submesh_1906-12-31T00:00:00"
nod=np.loadtxt(f"{M}/nod2d.out",skiprows=1); mp=np.loadtxt(f"{M}/map_nod.out").astype(int)
nl=np.loadtxt(f"{M}/cavity_nlvls.out").astype(int)
# zbar / Z levels from aux3d.out
aux=open(f"{M}/aux3d.out").read().split()
nlev=int(aux[0]); zbar=np.array(aux[1:1+nlev],float)
Z=0.5*(zbar[:-1]+zbar[1:])   # mid-depths
def rd(f,v):
    with nc.Dataset(f"{R}/{f}") as ds: return np.squeeze(np.asarray(ds.variables[v][:]))
T=rd("temp.nc","temp"); S=rd("salt.nc","salt")
def bymaxid(mid):
    return int(np.where(mp==mid)[0][0])
# crash element nodes: open crash 19690, deepened 19687, 19689
for mid in [19690,19687,19689,19668]:
    i=bymaxid(mid)
    print(f"\nnode maxid={mid} ulevels={nl[i]} lon{nod[i,1]:.2f} lat{nod[i,2]:.2f}")
    print("  lvl  Z(m)   T        S")
    for k in range(0,20):
        if T.ndim==2:
            print(f"   {k+1:2d} {Z[k]:6.0f} {T[k,i]:8.3f} {S[k,i]:8.3f}")
# horizontal T jump at matching levels between open crash (19690) and deepened (19689)
io=bymaxid(19690); ic=bymaxid(19689)
print("\n=== horizontal T/S jump: open 19690 vs deepened-cavity 19689 (matching level) ===")
print(" lvl  Z     T_open   T_cav   dT      S_open  S_cav   dS")
for k in range(6,18):
    print(f"  {k+1:2d} {Z[k]:5.0f} {T[k,io]:7.2f} {T[k,ic]:7.2f} {T[k,io]-T[k,ic]:7.2f}   {S[k,io]:6.2f} {S[k,ic]:6.2f} {S[k,io]-S[k,ic]:6.2f}")
# crude density (linear EOS) and column-integrated PGF proxy between the two columns
alphaT=0.2; betaS=0.8  # kg/m3 per degC / per psu (approx)
rho_o=1027 - alphaT*(T[:,io]-0) + betaS*(S[:,io]-34)
rho_c=1027 - alphaT*(T[:,ic]-0) + betaS*(S[:,ic]-34)
# only where both wet (level>=7)
g=9.81; L=12000.0
# integrate rho*g dz from top wet to bottom, difference / L / rho0 = accel scale
def colint(rho,ul):
    p=0.0; acc=[]
    for k in range(ul,18):
        p+=rho[k]*g*(zbar[k]-zbar[k+1])  # thickness>0
        acc.append(p)
    return np.array(acc)
ulc=nl[ic]-1
po=colint(rho_o,ulc); pc=colint(rho_c,ulc)
dP=po-pc
print(f"\nmax |dP| over shared column = {np.max(np.abs(dP)):.1f} Pa ; PGF accel ~ {np.max(np.abs(dP))/1027/L:.4e} m/s^2 ; *1200s = {np.max(np.abs(dP))/1027/L*1200:.3f} m/s")
