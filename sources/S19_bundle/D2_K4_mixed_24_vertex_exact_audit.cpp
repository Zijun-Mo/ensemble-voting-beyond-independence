// Certified complete vertex list of randomized/self-dual K=4 quadratic dual.
// Shows EVERY extreme ray of each of six vertex-orbit tangent cones either
// hits another one of 24 listed vertices or is an unbounded recession ray.
// Thus via connectedness of polyhedral 1-skeleton, no other vertices exist.
// Compile g++ -O3 -std=c++17 D2_K4_mixed_24_vertex_exact_audit.cpp -o audit
#include <bits/stdc++.h>
using namespace std;
using i64=long long; using i128=__int128_t;
constexpr int N=11,M=30;
using Z=array<i64,N>;
array<Z,M> A;array<i64,M> h;
vector<Z> V;vector<Z> P;
const int ee[6][2]={{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
int edge(int x,int y){if(x>y)swap(x,y);for(int e=0;e<6;e++)if(ee[e][0]==x&&ee[e][1]==y)return e;abort();}
Z permute(const Z&z,const array<int,4>&p){Z w{};w[0]=z[0];for(int i=0;i<4;i++)w[1+p[i]]=z[1+i];for(int k=0;k<6;k++)w[5+edge(p[ee[k][0]],p[ee[k][1]])]=z[5+k];return w;}
void setup(){
 int m=0;
 for(int mask=0;mask<16;mask++){
  Z row{};row[0]=1;for(int i=0;i<4;i++)row[1+i]=(mask>>i)&1;
  for(int k=0;k<6;k++)row[5+k]=((mask>>ee[k][0])&1)*((mask>>ee[k][1])&1);
  A[m++]=row;
 }
 for(int mask=0;mask<16;mask++)if(mask<(15^mask)){
  Z row{};for(int j=0;j<N;j++)row[j]=A[mask][j]+A[15^mask][j];
  A[m]=row;h[m++]=1;
 }
 for(int k=0;k<6;k++){Z row{};row[5+k]=1;A[m++]=row;}
 assert(m==M);
 P={Z{0,0,0,0,0,0,0,0,1,1,1},
    Z{0,0,0,0,1,0,0,0,0,0,0},
    Z{1,-1,-1,-1,-1,1,1,2,1,2,2},
    Z{1,-1,0,0,0,0,0,0,0,0,0},
    Z{3,-2,-2,-2,0,1,1,0,1,0,0},
    Z{6,-5,-3,-3,-3,2,2,2,1,1,1}};
 set<Z> uniq;
 array<int,4> p={0,1,2,3};
 for(const auto& prototype:P){do{uniq.insert(permute(prototype,p));}while(next_permutation(p.begin(),p.end()));p={0,1,2,3};}
 V.assign(uniq.begin(),uniq.end());
 assert(V.size()==24);
 for(auto const& v:V)for(int i=0;i<M;i++){
  i128 dot=0;for(int j=0;j<N;j++)dot+=(i128)A[i][j]*v[j];assert(dot>=h[i]);}
}
// A 10x11 full-row-rank integer matrix has a one-dimensional exact nullspace.
bool kernel(const array<int,10>&ii,Z &d){
 i64 T[10][11];for(int r=0;r<10;r++)for(int c=0;c<11;c++)T[r][c]=A[ii[r]][c];
 int pc[10];int rank=0; i64 prev=1;
 for(int c=0;c<11&&rank<10;c++){
  int s=rank;while(s<10&&T[s][c]==0)s++;
  if(s==10)continue;
  if(s!=rank)for(int t=c;t<11;t++)swap(T[rank][t],T[s][t]);
  i64 pivot=T[rank][c];
  for(int r=rank+1;r<10;r++){
   for(int t=c+1;t<11;t++){
    i128 num=(i128)T[r][t]*pivot-(i128)T[r][c]*T[rank][t];
    assert(num%prev==0);num/=prev;
    assert(num>=LLONG_MIN&&num<=LLONG_MAX);T[r][t]=(i64)num;
   }
   T[r][c]=0;
  }
  pc[rank]=c;prev=pivot;rank++;
 }
 if(rank<10)return false;
 int freecol=0;while(find(pc,pc+10,freecol)!=pc+10)freecol++;
 d.fill(0);d[freecol]=T[9][pc[9]];
 for(int r=9;r>=0;r--){
  i128 num=0;for(int c=pc[r]+1;c<11;c++)num-=(i128)T[r][c]*d[c];
  assert(num%T[r][pc[r]]==0);num/=T[r][pc[r]];
  assert(num>=LLONG_MIN&&num<=LLONG_MAX);d[pc[r]]=(i64)num;
 }
 i64 g=0;for(auto v:d)g=gcd(g,abs(v));assert(g>0);
 for(auto &v:d)v/=g;
 return true;
}
int active[30],kact,choose[10];
Z current;
long long bases=0,rank10=0,valid=0,bounded=0,unbounded=0;
set<Z> seen_dirs;
void classify(){
 ++bases;Z d;array<int,10>idx;for(int j=0;j<10;j++)idx[j]=active[choose[j]];
 if(!kernel(idx,d))return;
 ++rank10;
 bool plus=true,minus=true;
 for(int k=0;k<kact;k++){
  i128 v=0;for(int j=0;j<N;j++)v+=(i128)A[active[k]][j]*d[j];
  if(v<0)plus=false;
  if(v>0)minus=false;
 }
 if(!plus&&!minus)return;
 assert(plus!=minus);
 if(minus)for(auto &v:d)v=-v;
 if(!seen_dirs.insert(d).second)return;
 ++valid;
 // Follow this feasible tangent ray up to the nearest inactive facet.
 i64 mn=1,md=0; // fraction mn/md; md=0 denotes infinity
 for(int i=0;i<M;i++){
  i128 dot=0,slack=-h[i];
  for(int j=0;j<N;j++){
   dot+=(i128)A[i][j]*d[j];slack+=(i128)A[i][j]*current[j];
  }
  assert(slack>=0);
  if(dot<0){
   assert(slack>0);
   i64 nu=(i64)slack, de=(i64)-dot;
   if(md==0||(i128)nu*md<(i128)mn*de){mn=nu;md=de;}
  }
 }
 if(md==0){++unbounded;for(int i=0;i<M;i++){
  i128 dot=0;for(int j=0;j<N;j++)dot+=(i128)A[i][j]*d[j];assert(dot>=0);
 }return;}
 ++bounded;
 bool found=false;
 for(const Z&v:V){
  bool equal=true;
  for(int j=0;j<N;j++)if((i128)v[j]*md!=(i128)current[j]*md+(i128)mn*d[j]){equal=false;break;}
  if(equal){found=true;break;}
 }
 if(!found){cerr<<"UNKNOWN ADJACENT VERTEX at "<<" basis "<<bases<<"\n";exit(2);}
}
void comb(int start,int depth){
 if(depth==10){classify();return;}
 for(int k=start;k<=kact-(10-depth);k++){
  choose[depth]=k;comb(k+1,depth+1);
 }
}
int main(){
 setup();cout<<"EXACT: 24 candidates, all satisfy 30 inequalities.\n";
 for(int orbit=0;orbit<6;orbit++){
  current=P[orbit];kact=0;
  for(int i=0;i<M;i++){
   i128 v=-h[i];for(int j=0;j<N;j++)v+=(i128)A[i][j]*current[j];
   if(v==0)active[kact++]=i;
  }
  bases=rank10=valid=bounded=unbounded=0;seen_dirs.clear();
  comb(0,0);
  cout<<"orbit "<<orbit<<" act="<<kact<<" bases="<<bases<<" rank10="<<rank10
      <<" rays="<<valid<<" finite="<<bounded<<" infinity="<<unbounded<<"\n"<<flush;
 }
 cout<<"THEOREM: no missing vertices; 24 is the complete mixed dual vertex atlas.\n";
}
