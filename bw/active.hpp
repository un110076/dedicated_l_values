#pragma once

#include "dag.hpp"

/// construction of dag by overloading for custom data type
struct active { 
  /// primal function value
  double value=0; 
  /// reference to dag node
  int node_ref=-1; 

  active()=default;
  ~active()=default;
  active(const active&)=default;
  active(active&&)=default;
  active& operator=(const active&)=default;
  active& operator=(active&&)=default;
  active(const double& v) : value(v) {}

  /// ensures access to adjoint inputs
  void register_input() {
    G.num_indeps++;
    node_ref=G.dependencies.size();
    G.dependencies.push_back(G.adjoints_counter++);
  }

  /// ensures uniqueness of adjoint outputs and
  /// access to read-only inputs
  void register_output() {
    G.num_deps++;;
    // if (node_ref>=G.num_indeps) 
    //  G.dependencies[node_ref]=G.persistent_adjoints_counter--;
  }

  /// records index of argument of an operation and
  /// partial derivative of its result wrt. argument and
  /// updates bandwidth
  void record_arg(double d) const {
    G.bandwidth=std::max(G.bandwidth,G.adjoints_counter-G.dependencies[node_ref]);
    G.dependencies.push_back(G.dependencies[node_ref]);
    G.derivatives.push_back(d);
  }

  /// records number of arguments of an operation and
  /// index of its result
  void record_res(int num_args) {
    G.dependencies.push_back(num_args);
    node_ref=G.dependencies.size();
    G.dependencies.push_back(G.adjoints_counter++);
  }

};

#include "../active_ops.hpp"
