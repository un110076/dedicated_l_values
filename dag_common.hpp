#pragma once

  void dag::interpret() const {
    int node_ref=dependencies.size()-1;
    auto derivatives_it=derivatives.rbegin();
    while (node_ref>=num_indeps) {
      int res=node_ref--;
      int num_args=dependencies[node_ref--];
      double res_adj=adjoints[adjoint_idx(res)]; 
      adjoints[adjoint_idx(res)]=0;
      for (int i=0;i<num_args;++i)
        adjoints[adjoint_idx(node_ref--)]+=*derivatives_it++*res_adj;
      // print_adjoints();
    }
  }

  void dag::todot(const std::string& filename) const {
    using namespace std;
    ofstream out(filename);
    out << "digraph {" << endl << "rankdir=LR" << endl;
    int i=dependencies.size()-1;
    while (i>=0) {
      out << dependencies[--i] << ";" << endl;
      if (i>num_indeps) i-=dependencies[i]+1; 
    }
    i=dependencies.size()-1;
    while (i>=num_indeps) {
      int j=1;
      for (;j<=dependencies[i-1];j++) 
	out << dependencies[i-1-j] << "->" << dependencies[i] << ";" << endl;
      i-=j+1;
    }
    out << "}" << endl;
  }

  void dag::print() const {
    std::cerr << "dependencies" << std::endl;
    for (const auto& e:dependencies) 
      std::cerr << e << " ";
    std::cerr << std::endl;
    std::cerr << "derivatives" << std::endl;
    for (const auto& e:derivatives) 
      std::cerr << e << " ";
    std::cerr << std::endl;
  }
