import json, math, subprocess, os, sys
CA,CF = 3.0, 4.0/3.0
PY_ = "/usr/bin/python3"
ENV = dict(os.environ, PYTHONPATH="/home/cbiello/Packages/recola-heft-py/lib/python3.13/site-packages")
HERE = os.path.dirname(os.path.abspath(__file__))

def amp(**kw):
    r = subprocess.run([PY_, os.path.join(HERE,"amp.py"), json.dumps(kw)],
                       capture_output=True, text=True, env=ENV)
    try:   return json.loads(r.stdout.strip().split("\n")[-1])
    except Exception: print(r.stdout[-2000:], r.stderr[-2000:]); raise

def dot(a,b): return a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3]
def _boost(sign,q,p):
    """sign=+1: p seen from the rest frame of q;  sign=-1: inverse."""
    M=math.sqrt(dot(q,q)); g=q[0]/M
    b=[sign*q[k]/q[0] for k in (1,2,3)]
    bp=sum(b[k]*p[k+1] for k in range(3))
    out=[g*(p[0]-bp)]
    for k in range(3): out.append(p[k+1]+g*b[k]*(g/(g+1)*bp-p[0]))
    return out
def boost_to_rest(q,p):   return _boost(+1,q,p)
def boost_from_rest(q,p): return _boost(-1,q,p)

# ColombaPS massless point, POWHEG order 1,2 -> H b b~
P0=[[500.,0.,0.,500.],[500.,0.,0.,-500.],
    [211.07738134101575156,-144.92833555604772755,-30.63965724497451859,83.58020023094401552],
    [349.14239937112841972,-17.57847102882445967,-333.45082223855939674,101.99000707591839898],
    [439.78021928785585715,162.50680658487215169,364.09047948353401125,-185.57020730686238608]]

def massive_point(mb):
    """give legs 4,5 a mass mb keeping p4+p5, directions in its rest frame"""
    Q=[P0[3][i]+P0[4][i] for i in range(4)]
    s=dot(Q,Q); beta=math.sqrt(1-4*mb*mb/s); E=math.sqrt(s)/2
    out=[list(x) for x in P0]
    for leg in (3,4):
        q=boost_to_rest(Q,P0[leg])
        n=math.sqrt(sum(q[k]**2 for k in (1,2,3)))
        r=[E]+[q[k]/n*E*beta for k in (1,2,3)]
        out[leg]=boost_from_rest(Q,r)
    return out

def F1(mb,mu):                    # Mitov-Moch, two massive quark legs, as/2pi units
    L=math.log(mu*mu/(mb*mb))
    return CF*(L*L/2 + L/2 + 2 + math.pi**2/12)

mu, aS = 1000.0, 0.118
proc = sys.argv[1] if len(sys.argv)>1 else 'g g -> H b b~'
nf_massless = int(sys.argv[2]) if len(sys.argv)>2 else 5
nf_massive  = int(sys.argv[3]) if len(sys.argv)>3 else 4

m0 = amp(mb=0.0, nf=nf_massless, mu=mu, aS=aS, proc=proc, p=P0)
print("%s   mu=%g  nf(massless)=%d  nf(massive)=%d"%(proc,mu,nf_massless,nf_massive))
print("massless: f=%.8f  c1=%.8f  c2=%.8f"%(m0["f"],m0["p1"]-m0["f"],m0["p2"]-m0["f"]))
print()
print(" %8s %16s %16s %16s %14s"%("m_b","V/B massive","diff = m - 0","2*F1(MM)","diff - 2F1"))
for mb in (4.92,2.0,1.0,0.5,0.2,0.1,0.05,0.02):
    mm = amp(mb=mb, nf=nf_massive, mu=mu, aS=aS, proc=proc, p=massive_point(mb))
    d  = mm["f"]-m0["f"]
    print(" %8.3f %16.8f %16.8f %16.8f %14.8f"%(mb,mm["f"],d,F1(mb,mu),d-F1(mb,mu)))
