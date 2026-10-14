"""Self-contained rational certificate verifier for the heterogeneous K=4
closed form, using the finite vertex atlas and 262 Farkas certificates.

Usage: python D2_K4_exact_closed_form_verify.py
Optional: python D2_K4_exact_closed_form_verify.py --rebuild
The optional rebuild uses g++ to exhaust 705432 exact bases per orbit.
Requirements: Python standard library; g++ for optional --rebuild.
"""
import json,itertools,subprocess,sys,tempfile
from pathlib import Path
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parent
ATLAS=json.loads((ROOT/'D2_K4_all_signed_vertex_atlas.json').read_text())
WITNESS=json.loads((ROOT/'D2_K4_12term_Farkas_exact.json').read_text())
E=list(itertools.combinations(range(4),2))
X=[[int((m>>i)&1) for i in range(4)] for m in range(16)]
D=23

# Exact dual constraints: Q_z(x)>=F(x), gamma_ij >=0.
A=[[1,*x,*(x[i]*x[j] for i,j in E)] for x in X]
A += [[0]*5+[int(k==j) for k in range(6)] for j in range(6)]
C=[]
for i in range(4):
 c=[0]*11;c[1+i]=1;C.append(c)
for T in itertools.combinations(range(4),3):
 c=[0]*11
 for e in itertools.combinations(T,2):c[5+E.index(e)]=1
 C.append(c)
for r in range(4):
 c=[1]+[-1]*4+[2 if r in e else 1 for e in E]
 C.append(c)
assert len(C)==12 and len(A)==22

# Param domain x=(pi[16], b[6], t): pi>=0, b>=0,
# q=X^T pi <= 1/2, E[XiXj] <= b_ij <= min(q_i,q_j).
G=[];h=[]
for k,(i,j) in enumerate(E):
 for v in (i,j):
  z=[0]*D
  for m in range(16):z[m]=-X[m][v]
  z[16+k]=1;G.append(z);h.append(F(0))
 z=[0]*D
 for m in range(16):z[m]=X[m][i]*X[m][j]
 z[16+k]=-1;G.append(z);h.append(F(0))
for i in range(4):
 z=[0]*D
 for m in range(16):z[m]=X[m][i]
 G.append(z);h.append(F(1,2))
assert len(G)==22

def rational(s):return F(s)

def audit_record(rec):
 w=tuple(map(int,rec['w'].split(',')))
 den,*nn=rec['v']; v=[F(n,den) for n in nn]
 f=[int(sum(w[i]*(2*x[i]-1) for i in range(4))>0) for x in X]+[0]*6
 assert all(sum(F(a)*z for a,z in zip(row,v))>=rhs for row,rhs in zip(A,f))
 assert den>0
 # Add t<=L_C-L_V constraints to parameter polytope.
 rows=[list(map(F,z)) for z in G];hs=h.copy()
 for c in C:
  d=[F(u)-zz for u,zz in zip(c,v)]
  z=[F(0)]*D
  for m in range(16):z[m]=-sum(X[m][i]*d[1+i] for i in range(4))
  for k in range(6):z[16+k]=-d[5+k]
  z[22]=1
  rows.append(z);hs.append(d[0])
 assert len(rows)==34
 lamb=[F(0)]*34
 for i,s in rec['lam']:
  assert lamb[i]==0
  lamb[i]=rational(s)
 mu=rational(rec['mu'])
 assert all(z>=0 for z in lamb)
 # Stationarity for min (-t); pi,b>=0, t free.
 reduced=[F(int(j==22)*-1) for j in range(D)]
 for i,z in enumerate(lamb):
  if z:
   for j in range(D):reduced[j]+=z*rows[i][j]
 for j in range(16):reduced[j]+=mu
 assert all(z>=0 for z in reduced[:22])
 assert reduced[22]==0
 bound=-sum(a*b for a,b in zip(lamb,hs))-mu
 assert bound==rational(rec['bound'])>=0
 return (rec['w'],tuple(rec['v']))

seen=set()
for rec in WITNESS:
 key=audit_record(rec)
 assert key not in seen
 seen.add(key)
all_keys={(w,tuple(v)) for w,vertices in ATLAS.items() for v in vertices}
assert seen==all_keys
assert len(seen)==262 and len(ATLAS)==14
print('PASS: 262/262 exact rational Farkas inequalities; all 14 rule orbits')

# Independent exact basis enumeration, if requested.
if '--rebuild' in sys.argv:
 src=ROOT/'D2_K4_all_signed_rule_vertices_exact.cpp'
 with tempfile.TemporaryDirectory() as tmp:
  exe=Path(tmp)/'enumerate'
  subprocess.run(['g++','-O3','-std=c++17',str(src),'-o',str(exe)],check=True)
  for w,verts in ATLAS.items():
   p=subprocess.run([str(exe),*w.split(',')],check=True,capture_output=True,text=True)
   lines=p.stdout.splitlines()
   assert 'total_bases 705432 ' in lines[0]
   assert set(tuple(map(int,z.split())) for z in lines[1:])==set(tuple(v) for v in verts)
 print('PASS: independently rebuilt all 14 exact vertex atlases (14*705432 bases)')
print('THEOREM: R*_det,K4=min(q_i, triangle b-sums, H_r), all feasible q_i<1/2')
