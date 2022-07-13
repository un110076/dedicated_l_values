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
  int bandwidth=0, num_indeps=0, num_deps=0;
  static int adjoints_counter; 

  int RAM_size() {
    bandwidth=std::max(bandwidth,num_indeps);
    bandwidth=std::max(bandwidth,num_deps);
    return bandwidth;
  }

  int adjoint_idx(int dag_ref) const {
    int i=dependencies[dag_ref];
    if (i<0) return -i-1;
    return i%bandwidth;
  }

  void interpret() const;
  void todot(const std::string& filename) const;
  void print() const;

} G;

int dag::adjoints_counter=0;

#include "../dag_common.hpp"
