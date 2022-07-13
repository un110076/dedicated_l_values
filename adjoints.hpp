#pragma once

#include<vector>
#include<iostream>

static std::vector<double> adjoints;

static void print_adjoints() {
  std::cerr << "adjoints" << std::endl;
  for (const auto& e:adjoints)
    std::cerr << e << " ";
  std::cerr << std::endl;
}

