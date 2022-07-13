#include "active.hpp"

#include<vector>
#include<cassert>
#include<cmath>
#include<iostream>

template<typename T>
void f(std::vector<T>& v) {
  T u;
  for (size_t i=1;i<v.size();i++) { u=sin(v[i-1]); v[i]=u*u+v[0]; }
}

void fd_check(int l, double x) {
  using namespace std;
  vector<double> v(l);
  double h= (x==0)
    ? sqrt(numeric_limits<double>::epsilon())
    : sqrt(numeric_limits<double>::epsilon())*fabs(x);
  v[0]=x+h; f(v); double yp=v[l-1]; 
  v[0]=x-h; f(v); double ym=v[l-1];
  cout << "x_a=" << (yp-ym)/(2*h) << endl;
}

int main(int argc, char* argv[]) {
  using namespace std;
  assert(argc==3); 
  G.num_indeps=1;
  int l=stoi(argv[1]); assert(l>0);
  vector<active> v(l);
  double x=stof(argv[2]); 
  v[0]=x;
  v[0].register_input(); 
  int x_node_ref=v[0].node_ref;
  f(v);
  v[l-1].register_output(); 
  cout << "y=" << v[l-1].value << endl;
  G.todot("dag.dot");
  adjoints=vector<double>(G.RAM_size(),0);
  adjoints[G.adjoint_idx(v[l-1].node_ref)]=1;
  G.print(); 
  print_adjoints();
  G.interpret();
  print_adjoints();
  cerr << "size of RAM=" << adjoints.size()*8 << "b" << endl;
  cerr << "size of SAM=" << G.dependencies.size()*4+G.derivatives.size()*8 << "b" << endl;
  cout << "x_a=" << adjoints[G.adjoint_idx(x_node_ref)] << endl;
  cout << "CHECK:" << endl;
  fd_check(l,x);
  return 0;
}
