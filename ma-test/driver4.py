import json, math, subprocess, os, sys
D=os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(D,"driver.py")).read().split("mu, aS = ")[0])
mu, aS = 1000.0, 0.118
print("Loop content IDENTICAL on both sides (b massive with mass m in loops).")
print("massive side : g g -> H b b~ , external b mass m")
print("massless side: g g -> H s s~ , external s massless, MB = m in loops\n")
for proc_m, proc_0 in (('g g -> H b b~','g g -> H s s~'),
                       ('u u~ -> H b b~','u u~ -> H s s~')):
    print("%s   vs   %s"%(proc_m,proc_0))
    print("  %7s %9s %15s %15s %13s"%("m","L","diff","2*F1(MM)","diff-2*F1"))
    for mb in (4.92,1.0,0.2,0.05,0.02,0.01):
        mm = amp(mb=mb, nf=4, mu=mu, aS=aS, proc=proc_m, p=massive_point(mb))
        m0 = amp(mb=0.0, nf=4, mu=mu, aS=aS, proc=proc_0, p=P0, MBloop=mb)
        L=math.log(mu*mu/(mb*mb)); d=mm["f"]-m0["f"]
        print("  %7.3f %9.5f %15.8f %15.8f %13.8f"%(mb,L,d,2*F1(mb,mu),d-2*F1(mb,mu)))
    print()
