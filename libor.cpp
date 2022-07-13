#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <vector>
#include <random>
#include "active.hpp"

template<typename T>
void path_calc(
    const int N,
    const int Nmat,
    const T& delta,
    std::vector<T>& L,
    const std::vector<T>& sigma,
    const std::vector<double>& Z
) {
  using namespace std;
  for(int n=0;n<Nmat;n++) {
    T aux1=sqrt(delta)*Z[n];
    T S=0.0;
    for (int i=n+1;i<N;i++) {
      T aux2=delta*sigma[i-n-1];
      S=S+(aux2*L[i])/(1.0+delta*L[i]);
      L[i]=L[i]*exp(aux2*S+sigma[i-n-1]*(aux1-0.5*aux2));
    }
  }
}

template<typename T>
void portfolio(
    const int N,
    const int Nmat,
    const T& delta,
    const int Nopt,
    const std::vector<int>& maturities,
    const std::vector<T>& swaprates,
    const std::vector<T>& L,
    T& P
) {
  using namespace std;
  vector<T> B(N,0),S(N,0);
  T b=1.0;
  T s=0.0;
  for (int n=Nmat;n<N;n++) {
    b=b/(1.0+delta*L[n]); B[n]=b;
    s=s+delta*b; S[n]=s;
  }
  P=0;
  for (int i=0;i<Nopt;i++){
    int m=maturities[i]+Nmat-1;
    T swapval=B[m]+swaprates[i]*S[m]-1.0;
    if (swapval<0) P=P-100.0*swapval;
    // P=P+100.0*max(-swapval,0.0);
  }
  for (int n=0;n<Nmat;n++) P=P/(1.0+delta*L[n]);
}

template<typename T>
void libor(
    const int Nmat,
    const int Npath,
    const T& delta,
    const std::vector<int>& maturities,
    const std::vector<T>& swaprates,
    const std::vector<T>& sigma,
    const std::vector<T>& Lin,
    T& P
) {
  using namespace std;
  int N=Lin.size();
  int Nopt=maturities.size();
  vector<double> Z(Nmat,0);
  // vector<T> L(N,0);
  T Ps=0;
  // srand(0);
  default_random_engine generator(0);
  normal_distribution<double> distribution(0.0,1.0);

  for (int path=0; path<Npath; path++) {
    for (int i=0;i<Nmat;i++)
      Z[i]=0.3+distribution(generator);
    vector<T> L(Lin);
    path_calc(N,Nmat,delta,L,sigma,Z);
    portfolio(N,Nmat,delta,Nopt,maturities,swaprates,L,P);
    Ps=Ps+P;
  }
  P=Ps/Npath;
}

void fd_check(
    const int Nmat,
    const int Npath,
    const active delta,
    const std::vector<int>& maturities,
    const std::vector<active>& swaprates,
    const std::vector<active>& sigma,
    std::vector<active> Lin
) {
  using namespace std;
  int n=Lin.size();
  for (int i=0;i<n;i++) {
    active Pp=0,Pm=0;
    double h= (Lin[i].value==0)
      ? sqrt(numeric_limits<double>::epsilon())
      : sqrt(numeric_limits<double>::epsilon())*fabs(Lin[i].value);
    Lin[i].value+=h; 
    libor(Nmat,Npath,delta,maturities,swaprates,sigma,Lin,Pp);
    Lin[i].value-=2*h; 
    libor(Nmat,Npath,delta,maturities,swaprates,sigma,Lin,Pm);
    Lin[i].value+=h; 
    cout << "Lin_a[" << i << "]=" << (Pp.value-Pm.value)/(2*h) << endl;
  }
}

int main(int argc, char* argv[]) {
  using namespace std;
  assert(argc==2);
  int Npath=stoi(argv[1]); 
  const active delta=0.25;
  const vector<int> maturities({4,4,4,8,8,8,20,20,20,28,28,28,40,40,40});
  const vector<active> swaprates({.045,.05,.055,.045,.05,.055,.045,.05,.055,.045,.05,.055,.045,.05,.055 });
  const int Nmat=40;
  const int N=Nmat+40;
  vector<active> sigma(N,0.2);
  vector<active> L(N,0.05);
  active P=0;
  int n=N; G.num_indeps=n; 
  // register inputs
  vector<int> x_node_refs;
  for (auto& i:L) { i.register_input(); x_node_refs.push_back(i.node_ref); } 
  // record dag
  libor(Nmat,Npath,delta,maturities,swaprates,sigma,L,P);
  // register outputs
  P.register_output(); 
  // draw
  // G.todot("dag.dot");
  // allocate adjoints
  adjoints=vector<double>(G.RAM_size(),0);
  // seed
  adjoints[G.adjoint_idx(P.node_ref)]=1;
  // print
  // G.print(); 
  // print_adjoints();
  // interpret
  G.interpret();
  // print
  // print_adjoints();
  // results
  cerr << "size of RAM=" << adjoints.size()*8 << "b" << endl;
  cerr << "size of SAM=" << G.dependencies.size()*4+G.derivatives.size()*8 << "b" << endl;
  cout << "P=" << P.value << endl;
  for (int i=0;i<n;i++) 
    cout << "L_a[" << i << "]=" 
         << adjoints[G.adjoint_idx(x_node_refs[i])] << endl;
  // finite difference check
  // cout << "CHECK:"<< endl;
  // fd_check(Nmat,Npath,delta,maturities,swaprates,sigma,L);
  return 0;
}
