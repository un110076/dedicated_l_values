#pragma once

#include<vector>
#include<cassert>
#include<cmath>
#include<iostream>
#include<fstream>

#include "../adjoints.hpp"

static struct dag {
  std::vector<int> dependencies;
  std::vector<double> derivatives;
  int num_indeps, num_deps; 
  static int adjoints_counter; 

  int RAM_size() { return adjoints_counter; }

  int adjoint_idx(int node_ref) const {
    return dependencies[node_ref];
  }

  void interpret() const;
  void todot(const std::string& filename) const;
  void print() const;

} G;

int dag::adjoints_counter=0;

#include "../dag_common.hpp"
