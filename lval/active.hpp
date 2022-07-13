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
  active(const double& v) : value(v) {}

  /// ensures access to dedicated adjoint inputs and
  void register_input() {
    node_ref=G.dependencies.size();
    G.dependencies.push_back(G.persistent_adjoints_counter--);
  }

  /// dedicated outputs anyway
  void register_output() {}

  /// records index of argument of an operation and
  /// partial derivative of its result wrt. argument and
  /// updates bandwidth
  void record_arg(double d) const {
    if (G.dependencies[node_ref]>=0)
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

  /// dedicated adjoints for copies and non-move lhs
  void dedicate_copy_or_lhs(const active& arg) {
    arg.record_arg(1.);
    int d=node_ref;
    record_res(1);
    if (d<0)
      G.dependencies[node_ref]=G.persistent_adjoints_counter--;
    else
      G.dependencies[node_ref]=G.dependencies[d];
    G.adjoints_counter--;
  }

  /// dedicated adjoints for active copies 
  active(const active& arg) {
    if (arg.node_ref>=0) dedicate_copy_or_lhs(arg);
    value=arg.value;	
  }

  /// dedicated adjoints for active copies 
  // active(active&&)=default;
  // /*
  active(active&& arg) {
    if (arg.node_ref>=0) dedicate_copy_or_lhs(arg);
    value=arg.value;	
  }
  // */

  /// dedicated adjoints for active lhs
  active& operator=(const active& arg) {
    if (&arg!=this) {
      if (arg.node_ref>=0) dedicate_copy_or_lhs(arg);
      else node_ref=-1;
      value=arg.value;	
    }
    return *this;
  }

  /// dedicated adjoints for active lhs
  // active& operator=(active&&)=default;
  // /*
  active& operator=(active&& arg) {
    if (&arg!=this) {
      if (arg.node_ref>=0) dedicate_copy_or_lhs(arg);
      else node_ref=-1;
      value=arg.value;	
    }
    return *this;
  }
  // */
  
};

#include "../active_ops.hpp"
