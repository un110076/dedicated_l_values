#include<vector>
#include<cassert>
#include<cmath>
#include<iostream>
#include<random>

#include "active.hpp"

template<typename T>
void f(int np, const std::vector<T>& X, std::vector<T>& Y, const std::vector<double>& dW) {
  using namespace std;
  T s0=X[0];
  T e=X[1];
  T r=X[2];
  T sigma=X[3];
  T p=0;
  for (int i=0;i<np;i++) {
    T s=s0*exp((r-0.5*sigma*sigma)+sigma*dW[i]);
    p=p+exp(0-r)*max(s-e,0.0);
  }
  p=p/np;
  Y[0]=p;
}

void fd_check(void (*f)(int, const std::vector<double>&, std::vector<double>&, const std::vector<double>&), int n, int m, int np, const std::vector<double>& dW) {
  using namespace std;
  vector<double> x(n,2.), yp(m), ym(m); 
  for (int i=0;i<n;i++) {
    double h= (x[i]==0)
      ? sqrt(numeric_limits<double>::epsilon())
      : sqrt(numeric_limits<double>::epsilon())*fabs(x[i]);
    x[i]+=h; f(np,x,yp,dW); x[i]-=2*h; f(np,x,ym,dW);
    double r=0;
    for (int j=0;j<m;j++) r+=(yp[j]-ym[j])/(2*h);
    cout << "x_a[" << i << "]=" << r << endl;
  }
}

int main(int argc, char* argv[]) {
  using namespace std;
  assert(argc==2);
  int np=stoi(argv[1]);
  int n=4,m=1; G.num_indeps=n; 
  vector<active> x(n,2.), y(m); 
  srand(0);
  default_random_engine generator(0);
  normal_distribution<double> distribution(0.0,1.0);
  vector<double> dW(np);
  for (auto &dW_i:dW) dW_i=distribution(generator);
  // register inputs
  vector<int> x_node_refs;
  for (auto& i:x) { i.register_input(); x_node_refs.push_back(i.node_ref); } 
  // record dag
  f(np,x,y,dW);
  // register outputs
  for (auto& i:y) i.register_output();
  // draw
  // G.todot("dag.dot");
  // allocate adjoints
  adjoints=vector<double>(G.RAM_size(),0);
  // seed
  for (auto& i:y) adjoints[G.adjoint_idx(i.node_ref)]=1;
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
  for (size_t j=0;j<y.size();j++) 
    cout << "y[" << j << "]=" << y[j].value << endl;
  for (size_t i=0;i<x.size();i++) 
    cout << "x_a[" << i << "]=" << adjoints[G.adjoint_idx(x_node_refs[i])] << endl;
  // finite difference check
  //cout << "CHECK:"<< endl;
  //fd_check(f,n,m,np,dW);
  return 0;
}
