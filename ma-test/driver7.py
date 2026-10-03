import json, math, subprocess, os, sys
D=os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(D,"driver.py")).read().split("mu, aS = ")[0])
aS=0.118; mu=1000.; mb=0.02
def resid(Pm,P0m,proc_m,proc_0):
    Q=[Pm[3][i]+Pm[4][i] for i in range(4)]
    s=dot(Q,Q); beta=math.sqrt(1-4*mb*mb/s); E=math.sqrt(s)/2
    pm=[list(x) for x in Pm]
    for leg in (3,4):
        q=boost_to_rest(Q,Pm[leg]); n=math.sqrt(sum(q[k]**2 for k in (1,2,3)))
        pm[leg]=boost_from_rest(Q,[E]+[q[k]/n*E*beta for k in (1,2,3)])
    mm=amp(mb=mb,nf=4,mu=mu,aS=aS,proc=proc_m,p=pm)
    m0=amp(mb=0.0,nf=4,mu=mu,aS=aS,proc=proc_0,p=P0m,MBloop=mb)
    return mm["f"]-m0["f"]-2*F1(mb,mu)
Pswap=[P0[0],P0[1],P0[2],P0[4],P0[3]]      # b <-> bbar : same shat, same masses
print("same shat=1e6, same mu=1000, same m_b; only the SHAPE changes")
for pm,p0,tag in ((P0,P0,"original"),(Pswap,Pswap,"b <-> bbar")):
    for proc_m,proc_0 in (('g g -> H b b~','g g -> H s s~'),('u u~ -> H b b~','u u~ -> H s s~')):
        print("  %-12s %-16s residual = %14.8f"%(tag,proc_m,resid(pm,p0,proc_m,proc_0)))
