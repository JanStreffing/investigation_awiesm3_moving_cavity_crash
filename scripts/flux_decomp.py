import numpy as np, netCDF4 as nc
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
W="/work/ab0246/a270092/pism_repro/experiments/movcav4/run_awiesm3_19010101-19011231/work"
g=nc.Dataset(f"{W}/grids.nc"); alon=np.array(g["A096.lon"]).ravel(); alat=np.array(g["A096.lat"]).ravel(); g.close()
HLON,HLAT=108.73,-66.09
ia=int(np.argmin((np.cos(np.radians(HLAT))*(alon-HLON))**2+(alat-HLAT)**2))
def rd(f,v): d=nc.Dataset(f"{W}/{f}"); a=np.array(d[v])[:,0,:]; d.close(); return a
Qs =rd("A_Qs_all_OpenIFS_01.nc","A_Qs_all")    # solar to ocean
Qns=rd("A_Qns_oce_OpenIFS_01.nc","A_Qns_oce")  # non-solar (LW+sens+lat) to ocean
Qi =rd("A_Q_ice_OpenIFS_03.nc","A_Q_ice")      # to ice
ba=(alat<-64)&(alat>-68)&(alon>100)&(alon<118)
# coastal-max cell of NON-SOLAR (the suspicious one) at final step
tf=Qns.shape[0]-1
icmax=np.where(ba)[0][np.argmax(np.abs(Qns[tf,ba]))]
print(f"crash A096 idx={ia} @({alon[ia]:.2f},{alat[ia]:.2f})")
print(f"coastal max-|Qns| cell idx={icmax} @({alon[icmax]:.2f},{alat[icmax]:.2f})")
print(f"\n           |   solar Qs        |  non-solar Qns    |   to-ice Qi")
print(f"cell       |  mean   min   max |  mean   min   max |  mean   max")
for nm,i in [("crash",ia),("maxQns",icmax)]:
    print(f"{nm:10s} | {Qs[:,i].mean():5.0f} {Qs[:,i].min():5.0f} {Qs[:,i].max():5.0f} |"
          f" {Qns[:,i].mean():5.0f} {Qns[:,i].min():5.0f} {Qns[:,i].max():5.0f} |"
          f" {Qi[:,i].mean():5.0f} {Qi[:,i].max():5.0f}")
# regional peak per component
print(f"\nregional (coast) peak: Qs={np.nanmax(Qs[:,ba]):.0f}  |Qns|={np.nanmax(np.abs(Qns[:,ba])):.0f}  |Qi|={np.nanmax(np.abs(Qi[:,ba])):.0f} W/m2")
fig,ax=plt.subplots(2,1,figsize=(11,7),sharex=True)
for k,(nm,i) in enumerate([("crash cell",ia),("coastal max-|Qns| cell",icmax)]):
    ax[k].plot(Qs[:,i],'orange',lw=1.6,label='solar Qs (A_Qs_all)')
    ax[k].plot(Qns[:,i],'purple',lw=1.6,label='non-solar Qns (A_Qns_oce: LW+sens+lat)')
    ax[k].plot(Qs[:,i]+Qns[:,i],'k--',lw=1,alpha=.6,label='total to ocean')
    ax[k].axhline(0,color='gray',lw=.6); ax[k].set_ylabel('W/m2'); ax[k].set_title(f"{nm} @({alon[i]:.1f}E,{alat[i]:.1f}N)")
    ax[k].legend(fontsize=8,loc='upper right')
ax[1].set_xlabel('coupling step (hourly)')
plt.suptitle('OIFS->ocean heat flux decomposition: solar vs non-solar (Jan, SH summer)')
plt.tight_layout(); plt.savefig("/home/a/a270092/esm_tools/flux_decomp_solar_vs_nonsolar.png",dpi=140)
print("\nsaved flux_decomp_solar_vs_nonsolar.png")
