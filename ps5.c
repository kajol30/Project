#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define K 5
int N,M,*L,*occ,*os,*tc,*cr,*un,*up,nu;char*v;
unsigned long r=88172645463325252UL;static inline unsigned long xr(){r^=r<<13;r^=r>>7;r^=r<<17;return r;}
static inline int sat(int l){return l>0?v[l]:!v[-l];}
void addU(int c){up[c]=nu;un[nu++]=c;} void remU(int c){int p=up[c],l=un[--nu];un[p]=l;up[l]=p;}
int main(int a,char**g){double lim=atof(g[1]);r^=atol(g[2])*2654435761UL;r|=1;double cb=a>3?atof(g[3]):3.7;
FILE*f=fopen(a>4?g[4]:"cnf384.txt","r");fscanf(f,"p cnf %d %d",&N,&M);L=malloc(M*K*4);int*cnt=calloc(N+2,4);
for(int i=0;i<M;i++){for(int k=0;k<K;k++){fscanf(f,"%d",&L[i*K+k]);cnt[abs(L[i*K+k])]++;}int z;fscanf(f,"%d",&z);}
os=malloc((N+2)*4);int s=0;for(int x=1;x<=N+1;x++){os[x]=s;if(x<=N)s+=cnt[x];}occ=malloc(s*4);int*fi=calloc(N+2,4);
for(int i=0;i<M;i++)for(int k=0;k<K;k++){int x=abs(L[i*K+k]);occ[os[x]+fi[x]++]=i;}
v=malloc(N+1);tc=malloc(M*4);cr=malloc(M*4);un=malloc(M*4);up=malloc(M*4);
double PT[64];for(int i=0;i<64;i++)PT[i]=pow(cb,-i);
clock_t t0=clock();int best=M;
while((clock()-t0)/(double)CLOCKS_PER_SEC<lim){
 for(int x=1;x<=N;x++)v[x]=xr()&1;nu=0;
 for(int c=0;c<M;c++){tc[c]=0;cr[c]=0;for(int k=0;k<K;k++)if(sat(L[c*K+k])){tc[c]++;cr[c]^=abs(L[c*K+k]);}if(!tc[c])addU(c);}
 for(long st=0;st<20000000L&&nu;st++){
  int c=un[xr()%nu];double pr[K],sm=0;
  for(int k=0;k<K;k++){int x=abs(L[c*K+k]);int b=0;for(int j=os[x];j<os[x+1];j++){int d=occ[j];if(tc[d]==1&&cr[d]==x)b++;}pr[k]=PT[b<63?b:63];sm+=pr[k];}
  double rr=(xr()%1000000)/1e6*sm;int k=0;while(k<K-1&&rr>=pr[k]){rr-=pr[k];k++;}
  int x=abs(L[c*K+k]);v[x]^=1;
  for(int j=os[x];j<os[x+1];j++){int d=occ[j];int lt=0;for(int q=0;q<K;q++)if(abs(L[d*K+q])==x){lt=L[d*K+q];break;}
   if(sat(lt)){tc[d]++;cr[d]^=x;if(tc[d]==1)remU(d);}else{tc[d]--;cr[d]^=x;if(tc[d]==0)addU(d);}}
  if(nu<best){best=nu;}
  if((st&0xFFFF)==0&&(clock()-t0)/(double)CLOCKS_PER_SEC>lim)break;}
 if(!nu){for(int x=1;x<=N;x++)putchar('0'+v[x]);putchar('\n');return 0;}}
fprintf(stderr,"best %d\n",best);return 1;}
