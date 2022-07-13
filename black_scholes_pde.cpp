#include<vector>
#include<cassert>
#include<cmath>
#include<iostream>

#include "active.hpp"

template<typename T>
void f(std::vector<T> x, std::vector<T>& y) {
  const T& e=x[x.size()-3];
  const T& r=x[x.size()-2];
  const T& sigma=x[x.size()-1];
  std::vector<T>& u=x;
  std::vector<T> u_new(u.size()-3);
  int nx=u.size()-2, nt=nx*nx;
  double dt=1./nt;
  for (int j=0;j<nt;j++) {
    double t=j*dt;
    u_new[0]=u[0]+dt*(
      0.5*pow(sigma*(0+1),2)*(u[0+1]-2*u[0]+/*=u[0-1]=*/0/**/)+
      0.5*r*(0+1)*(u[0+1]-/*u[0-1]=*/0/**/)-r*u[0]);
    for (int i=1;i<nx-2;i++)
      u_new[i]=u[i]+dt*
        (0.5*pow(sigma*(i+1),2)*(u[i+1]-2*u[i]+u[i-1])+
        0.5*r*(i+1)*(u[i+1]-u[i-1])-r*u[i]);
    u_new[nx-2]=u[nx-2]+dt*(
      0.5*pow(sigma*(nx-2+1),2)*(/*u[nx-2+1]=*/1-e*exp(-r*t)/**/-2*u[nx-2]+u[nx-2-1])+
      0.5*r*(nx-2+1)*(/*u[nx-2+1]=*/1-e*exp(-r*t)/**/-u[nx-2-1])-r*u[nx-2]);
    u=u_new;
  }
  y=u_new;
}

void fd_check(void (*f)(std::vector<double>, std::vector<double>&), int n, int m) {
  using namespace std;
  vector<double> x(n,.2), yp(m), ym(m); 
  for (int i=0;i<n;i++) {
    double h= (x[i]==0)
      ? sqrt(numeric_limits<double>::epsilon())
      : sqrt(numeric_limits<double>::epsilon())*fabs(x[i]);
    x[i]+=h; f(x,yp); x[i]-=2*h; f(x,ym);
    double r=0;
    for (int j=0;j<m;j++) r+=(yp[j]-ym[j])/(2*h);
    cout << "x_a[" << i << "]=" << r << endl;
  }
}

int main(int argc, char* argv[]) {
  using namespace std;
  assert(argc==2);
  int n=stoi(argv[1]);
  int m=n-3; 
  G.num_indeps=n; 
  vector<active> x(n,.2), y(m); 
  // register inputs
  vector<int> x_node_refs;
  for (auto& i:x) { i.register_input(); x_node_refs.push_back(i.node_ref); }
  // record dag
  f(x,y);
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
    cout << "x_a[" << i << "]=" 
         << adjoints[G.adjoint_idx(x_node_refs[i])] << endl;
  // finite difference check
  // cout << "CHECK:"<< endl;
  // fd_check(f,n,m);
  return 0;
}
