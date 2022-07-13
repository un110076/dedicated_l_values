#include <cstdlib>
#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
using namespace std;

#include "active.hpp"

template <typename T>
inline T norm(vector<T>& v) {
  int n=v.size();
  T r=0;
  for (int i=1;i<n-1;i++) r=r+v[i]*v[i];
  return sqrt(r);
}

const double pi=3.141592653589793;

struct progress {
  static void update(float progress, bool output=true) {
    if (!output) return;
    int barWidth = 70;

    std::cerr << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
      if (i < pos) std::cerr << "=";
      else if (i == pos) std::cerr << ">";
      else std::cerr << " ";
    }
    std::cerr << "] " << int(progress * 100.0) << " %\r";
    std::cerr.flush();
  }
  static void finish(bool output=true) {
    if (!output) return;
    update(1.0);
    std::cerr << std::endl;
  }
};

// A:=L+U-I s.t. A=L*U
template <typename T>
inline void LU(vector<T>& A) {
  int n=(A.size()-4)/3+2;
  for (int k=1;k<n-1;k++) {
    A[2+k*3]=A[2+k*3]/A[k*3];
    A[(k+1)*3]=A[(k+1)*3]-A[2+k*3]*A[(k+1)*3-2];
  }
}

// y:=L^-1*b
template <typename T>
inline void FS(const vector<T>& LU, vector<T>& b) {
  int n=b.size();
  for (int i=1;i<n-1;i++)
    b[i]=b[i]-LU[i*3-1]*b[i-1];
}

// x:=U^-1*y
template <typename T>
inline void BS(const vector<T>& LU, vector<T>& y) {
  int n=y.size();
  for (int k=n-1;k>1;k--) {
    y[k-1]=y[k-1]-LU[(k-1)*3+1]*y[k];
    y[k-1]=y[k-1]/LU[(k-1)*3];
  }
}


//*** right-hand side of ODE dy/dt = v d^2 y / d x^2 - y dy / dx
//***                                `---,---------´   `---,---´
//***                                    diffusion         advection
template <typename T>
inline void g(const double& d, const vector<T>& y, vector<T>& r) {
  int n=y.size();
  for (int i=1;i<n-1;i++) {
    //*** central finite difference scheme for diffusion term on
    //*** equidistant 1D grid
    T diffusion = d*(n-1)*(n-1)*(y[i-1]-2*y[i]+y[i+1]);
    //*** first-order upwind scheme for advection term on equidistant
    //*** 1D grid
    T advection = -y[i]*(n-1);
    if (advection < 0) { advection = advection*(y[i]-y[i-1]); }
    else               { advection = advection*(y[i+1]-y[i]); }
    r[i] = diffusion + advection;
  }
}

// tangent of g (hand-written)
template <typename T>
inline void g_t(const double& d, const vector<T>& y, const vector<T>& y_t, vector<T>& r_t) {
  int n=y.size();
  for (int i=1;i<n-1;i++)  {
    T diffusion_t = d*(n-1)*(n-1)*(y_t[i-1]-2*y_t[i]+y_t[i+1]);
    T advection = -y[i]*(n-1);
    T advection_t = -y_t[i]*(n-1);
    if (advection < 0) { advection_t = advection * (y_t[i]-y_t[i-1]) + advection_t * (y[i]-y[i-1]); }
    else               { advection_t = advection * (y_t[i+1]-y_t[i]) + advection_t * (y[i+1]-y[i]); }
    r_t[i] = diffusion_t + advection_t;
  }
}

//*** Jacobian of r with respect to y (i.e. Jacobian of implementation g)
//***  - for the computation, the tangent of g (i.e. g_t) is used
//***  - A is known to be a tridiagonal matrix and uses custom data format, it stores row-wise
//***       a b
//***       c d e
//***         f g h
//***           i j k
//***             l m
//***    where in our case (a, b, c, k, l, m) = 0
//***  - Curtis/Powell/Reid(CPR)-algorithm is used to compute compressed Jacobian
//***  - if transpose=true, the transposed Jacobian is returned in A using the same data format.
template <typename T>
inline void dgdy(const double& d, const vector<T>& y, vector<T>& A,
                 bool transpose=false) {
  int n=y.size();
  for (int i=0;i<3;i++) {
    vector<T> y_t(n,0),r_t(n);
    //*** CPR seeding
    for (int j=i+1;j<n;j+=3) y_t[j]=1;
    
    //*** computes compressed column in r_t
    g_t(d,y,y_t,r_t);

    for (int j=i+1;j<n-1;j+=3) {
      A[j*3]=r_t[j]; // diagonal elements (e.g. g)
      if (transpose) {
        //*** store column r_t in row
        A[j*3-1]=r_t[j-1]; // left of diagonal (i.e. f)
        A[j*3+1]=r_t[j+1]; // right of diagonal (i.e. h)
      } else {
        A[(j-1)*3+1]=r_t[j-1]; // above diagonal (i.e. e)
        A[(j+1)*3-1]=r_t[j+1]; // below diagonal (i.e. i)
      }
    }
  }
}

//*** residual of nonlinear system
//*** - y_prev is solution of previous timestep
template <typename T>
inline void f(const int m, const double& d,
    const vector<T>& y, const vector<T>& y_prev, vector<T>& r) {
  int n=y.size();
  g(d,y,r);
  r[0]=y[0];
  for (int i=1;i<n-1;i++) r[i]=y[i]-y_prev[i]-r[i]/m;
  r[n-1]=y[n-1];
}

//*** Jacobian of residual of nls wrt. state (i.e. y)
//*** - see description of dgdy
template <typename T>
inline void dfdy(const int m, const double& d,
                 const vector<T>& y,
                 vector<T>& A, bool transpose=false) {
  int n=y.size();
  dgdy(d,y,A,transpose);
  for (auto& e : A) e=-e/m;
  for (int i=1;i<n-1;i++) A[i*3]=A[i*3]+1; 
}

//*** Newton solver for nls
template <typename T>
inline void newton(const int m, const double& d, const vector<T>& y_prev, vector<T>& y) {
  int n=y.size();
  const double eps=1e-14;
  //*** see comment of dgdy for description of data structure A
  vector<T> A((n-2)*3+4,0), r(n,0);
  f(m,d,y,y_prev,r); 
  while (norm(r)>eps) {
    dfdy(m,d,y,A);
    LU(A); FS(A,r); BS(A,r);
    for (int i=1;i<n-1;i++) y[i]=y[i]-r[i]; 
    f(m,d,y,y_prev,r); 
  }
}

//*** implicit Euler integration
template <typename T>
inline void burgers(const int m, const double& d, vector<T>& y, bool output_progress=true) {
  for (int j=0;j<m;j++) {
//    progress::update(static_cast<float>(j)/static_cast<float>(m), output_progress);
    vector<T> y_prev=y;
    newton(m,d,y_prev,y); 
  }
//  progress::finish(output_progress);
}

int main(int argc, char* argv[]){
  assert(argc==3);
  int n=stoi(argv[1]), m=stoi(argv[2]);
  G.num_indeps=n;
  vector<active> y(n,0);
  for (int i=1;i<n-1;i++) y[i]=sin((2*pi*i)/n);
  const double d=1e-3;
  // register inputs
  vector<int> y_node_refs;
  for (auto& i:y) { i.register_input(); y_node_refs.push_back(i.node_ref); }
  // record dag
  burgers(m,d,y);
  // register outputs
  for (auto& i:y) i.register_output();
  // draw
  // G.todot("dag.dot");
  // allocate adjoints
  adjoints=vector<double>(G.RAM_size(),0);
  // seed
  for (int i=1;i<n-1;i++) adjoints[G.adjoint_idx(y[i].node_ref)]=1; 
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
  for (size_t i=0;i<y.size();i++)
    cout << "y_a[" << i << "]=" 
         << adjoints[G.adjoint_idx(y_node_refs[i])] << endl;
  return 0;
}
