import json, math, subprocess, os, sys
D=os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(D,"driver.py")).read().split("mu, aS = ")[0])
aS=0.118
print("residual = [V/B]_massive - [V/B]_massless - 2*F1(MitovMoch)   (eps^0, Delta=0)")
print("%-16s %6s %7s  %14s"%("process","lambda","mu","residual"))
for proc_m,proc_0 in (('g g -> H b b~','g g -> H s s~'),('u u~ -> H b b~','u u~ -> H s s~')):
  for lam,mu in ((1.0,1000.),(3.0,1000.),(0.5,1000.),(1.0,300.),(1.0,3000.)):
    Pl=[[lam*x for x in r] for r in P0]; MH=125.*lam
    mb=0.02
    Q=[Pl[3][i]+Pl[4][i] for i in range(4)]
    s=dot(Q,Q); beta=math.sqrt(1-4*mb*mb/s); E=math.sqrt(s)/2
    pm=[list(x) for x in Pl]
    for leg in (3,4):
        q=boost_to_rest(Q,Pl[leg]); n=math.sqrt(sum(q[k]**2 for k in (1,2,3)))
        pm[leg]=boost_from_rest(Q,[E]+[q[k]/n*E*beta for k in (1,2,3)])
    mm=amp(mb=mb,nf=4,mu=mu,aS=aS,proc=proc_m,p=pm,MH=MH)
    m0=amp(mb=0.0,nf=4,mu=mu,aS=aS,proc=proc_0,p=Pl,MH=MH,MBloop=mb)
    print("%-16s %6.1f %7.0f  %14.8f"%(proc_m,lam,mu,mm["f"]-m0["f"]-2*F1(mb,mu)))
