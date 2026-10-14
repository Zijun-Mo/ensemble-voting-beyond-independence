// Exhaust all C(22,11)=705432 active bases of the K=4 anchor-rule LP dual.
// Fraction-free Bareiss elimination and exact integer feasibility tests.
// Compile: g++ -O3 -std=c++17 D2_K4_anchor_all_vertices_exact.cpp -o k4_vertices
// Run: ./k4_vertices
#include <array>
#include <set>
#include <vector>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <cassert>
using i64=long long;
using i128=__int128_t;
using Vec=std::array<i64,12>;
constexpr int N=11, M=22;
std::array<Vec,M> rows;
std::array<int,N> chosen;
std::set<std::array<i64,N+1>> vertices; // common positive denominator plus 11 numerators
long long all_bases=0,nonsingular_bases=0,feasible_bases=0;

void solve_basis(){
    ++all_bases;
    i64 T[N][N+1];
    for(int i=0;i<N;i++) for(int j=0;j<=N;j++) T[i][j]=rows[chosen[i]][j];
    i64 previous=1;
    for(int k=0;k<N;k++){
        int pivot=k;
        while(pivot<N && T[pivot][k]==0) ++pivot;
        if(pivot==N) return;
        if(pivot!=k) for(int j=k;j<=N;j++) std::swap(T[k][j],T[pivot][j]);
        i64 d=T[k][k];
        if(k<N-1){
            for(int i=k+1;i<N;i++){
                for(int j=k+1;j<=N;j++){
                    i128 num=(i128)T[i][j]*d-(i128)T[i][k]*T[k][j];
                    assert(num%previous==0);
                    i128 v=num/previous;
                    assert(v>=INT64_MIN && v<=INT64_MAX);
                    T[i][j]=(i64)v;
                }
                T[i][k]=0;
            }
        }
        previous=d;
    }
    ++nonsingular_bases;
    i64 denom=T[N-1][N-1];
    i64 numer[N];
    for(int j=N-1;j>=0;j--){
        i128 val=(i128)T[j][N]*denom;
        for(int k=j+1;k<N;k++) val-=(i128)T[j][k]*numer[k];
        assert(val%T[j][j]==0);
        val/=T[j][j];
        assert(val>=INT64_MIN && val<=INT64_MAX);
        numer[j]=(i64)val;
    }
    if(denom<0){denom=-denom;for(int j=0;j<N;j++)numer[j]=-numer[j];}
    // Full dual feasibility with exact arithmetic.
    for(int i=0;i<M;i++){
        i128 lhs=0;
        for(int j=0;j<N;j++)lhs+=(i128)rows[i][j]*numer[j];
        if(lhs<(i128)rows[i][N]*denom)return;
    }
    ++feasible_bases;
    i64 gcd=denom;
    for(int j=0;j<N;j++)gcd=std::gcd(gcd,numer[j]);
    std::array<i64,N+1> key{};
    key[0]=denom/gcd;
    for(int j=0;j<N;j++)key[j+1]=numer[j]/gcd;
    vertices.insert(key);
}
void combinations(int start,int pos){
    if(pos==N){solve_basis();return;}
    for(int i=start;i<=M-(N-pos);i++){
        chosen[pos]=i;combinations(i+1,pos+1);
    }
}
int main(int argc,char**argv){
    int w[4]={2,1,1,1};
    if(argc==5)for(int j=0;j<4;j++)w[j]=std::stoi(argv[j+1]);
    // Variables a,beta_0..beta_3,gamma_01,gamma_02,gamma_03,
    // gamma_12,gamma_13,gamma_23.
    constexpr int e1[6]={0,0,0,1,1,2};
    constexpr int e2[6]={1,2,3,2,3,3};
    for(int mask=0;mask<16;mask++){
        auto &r=rows[mask];r.fill(0);r[0]=1;
        int x[4];for(int i=0;i<4;i++){x[i]=(mask>>i)&1;r[1+i]=x[i];}
        for(int k=0;k<6;k++)r[5+k]=x[e1[k]]*x[e2[k]];
        int score=0;for(int i=0;i<4;i++)score+=w[i]*(2*x[i]-1);
        assert(score!=0);r[11]=(score>0)?1:0;
    }
    for(int k=0;k<6;k++){
        auto &r=rows[16+k];r.fill(0);r[5+k]=1;
    }
    combinations(0,0);
    std::cout<<"total_bases "<<all_bases<<" nonsingular "<<nonsingular_bases
             <<" feasible_bases "<<feasible_bases<<" distinct_vertices "<<vertices.size()<<"\n";
    assert(all_bases==705432);
    for(auto const &v:vertices){
        for(int i=0;i<N+1;i++) std::cout<<(i?" ":"")<<v[i];
        std::cout<<"\n";
    }
}
