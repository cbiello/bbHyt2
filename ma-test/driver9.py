import json, math, os
D=os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(D,"driver.py")).read().split("mu, aS = ")[0])
aS=0.118; z2=math.pi**2/6
def F1ttH(mb,mu):                   # ttH_NLO massif_factors, 2 legs, as/2pi
    L=math.log(mu*mu/(mb*mb)); return CF*(L*L/2+L/2+2+math.pi**2/12)
print("PRODUCTION setup: massless side is the SAME process with MB=0,")
print("so the b is a massless flavour in the loops (nf=5), while the")
print("exact run has it decoupled (Nf4).")
print("residual = V^m_BLHA/B - [2*F1_ttH + V^0_BLHA/B + C_F*pi^2/6]\n")
for proc in ('g g -> H b b~','u u~ -> H b b~'):
  for mu in (1000.,3000.):
    print("  %-16s mu=%5.0f"%(proc,mu))
    rows=[]
    for mb in (0.2,0.05,0.01):
        mm=amp(mb=mb,nf=4,mu=mu,aS=aS,proc=proc,p=massive_point(mb))
        m0=amp(mb=0.0,nf=5,mu=mu,aS=aS,proc=proc,p=P0)
        VB_m = mm["f"]+(mm["p2"]-mm["f"])*z2
        VB_0 = m0["f"]+(m0["p2"]-m0["f"])*z2
        L=math.log(mu*mu/(mb*mb))
        r = VB_m - (2*F1ttH(mb,mu) + VB_0 + CF*math.pi**2/6)
        rows.append((L,r)); print("     m=%-6g L=%8.5f  residual=%12.7f"%(mb,L,r))
    (L1,r1),(L2,r2)=rows[0],rows[-1]
    sl=(r2-r1)/(L2-L1); print("     slope=%.6f  const=%.6f"%(sl,r1-sl*L1))
