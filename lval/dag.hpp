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
  int bandwidth=0, num_indeps, num_deps;
  static int adjoints_counter; 
  static int persistent_adjoints_counter; 
  static int RAM_size_et; 
  static int RAM_size_noet; 

  int RAM_size() const {
    return bandwidth-persistent_adjoints_counter-1;
  }

  int adjoint_idx(int dag_ref) const {
    int i=dependencies[dag_ref];
    if (i<0) return -i-1;
    return i%bandwidth-persistent_adjoints_counter-1;
  }

  void interpret() const;
  void todot(const std::string& filename) const;
  void print() const;

} G;

int dag::adjoints_counter=0;
int dag::persistent_adjoints_counter=-1;
int dag::RAM_size_et=0;
int dag::RAM_size_noet=0;

#include "../dag_common.hpp"
