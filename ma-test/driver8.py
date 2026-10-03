import json, math, os
D=os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(D,"driver.py")).read().split("mu, aS = ")[0])
aS=0.118; mu=1000.; z2=math.pi**2/6
print("massification constant 2*F1 extracted in the two Recola conventions")
print("  Delta=0 : expect C_F(L^2+L+4)          = ttH form - C_F*pi^2/6")
print("  BLHA    : expect C_F(L^2+L+4+pi^2/3)   = ttH form + C_F*pi^2/6\n")
print("  %-16s %-8s %14s %14s"%("process","m_b","Delta=0","BLHA"))
for proc_m,proc_0 in (('g g -> H b b~','g g -> H s s~'),('u u~ -> H b b~','u u~ -> H s s~')):
  for mb in (0.05,0.01):
    mm=amp(mb=mb,nf=4,mu=mu,aS=aS,proc=proc_m,p=massive_point(mb))
    m0=amp(mb=0.0,nf=4,mu=mu,aS=aS,proc=proc_0,p=P0,MBloop=mb)
    L=math.log(mu*mu/(mb*mb))
    d0 = mm["f"]-m0["f"]
    dB = (mm["f"]+(mm["p2"]-mm["f"])*z2) - (m0["f"]+(m0["p2"]-m0["f"])*z2)
    print("  %-16s %-8g %14.8f %14.8f"%(proc_m,mb,d0-CF*(L*L+L+4),dB-CF*(L*L+L+4+math.pi**2/3)))
