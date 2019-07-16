#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include "fista.hpp"
#include <iostream>
#include <cmath>

namespace dtd {
  void invxtgx(mat const & x, vec const & g, mat & xi) {
    mat xtgx = x.transpose()*g.asDiagonal()*x;
    xi.resize(x.cols(), x.cols());
    xi = xtgx.llt().solve(mat::Identity(x.cols(), x.cols()));
    std::cout <<"from c++:\n"<< xi << "\n";
    // return xi;
  }
  namespace stat {
    ftype var(vec const & v) {
      vec vm = v.array() - v.mean();
      auto n = v.size();
      return vm.dot(vm) / (n-1);
    }
    ftype cov(vec const & a, vec const & b) {
      int n = a.size();
      assert(n == b.size());
      vec am = a.array() - a.mean();
      vec bm = b.array() - b.mean();
      return am.dot(bm)/(n-1);
    }
  }
  ftype GoertlerModel::eval(vec const & g) const {
    return 0.0;
  }
  void GoertlerModel::grad(vec & gr, vec const & g) const {
  }

  template<class Model>
  double FistaSolver<Model>::iterate() {
    std::cout << "hello from FistaSolver::iterate()\n";
    return 0.0;
  }
}
