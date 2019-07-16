#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include "fista.hpp"
#include <iostream>
#include <cmath>

namespace dtd {
  mat invxtgx(mat const & x, vec const & g) {
    mat xtgx = x.transpose()*g.asDiagonal()*x;
    mat xi = xtgx.llt().solve(mat::Identity(x.cols(), x.cols()));
    std::cout << xi << "\n";
    return xi;
  }
  namespace stat {
    ftype mean(vec const & v) {
      return v.sum() / v.size();
    }
    ftype var(vec const & v) {
      ftype f = mean(v);
      ftype res(0.0);
      for( int i = 0; i < v.size(); ++i) {
        res += std::pow((f - v(i)), 2);
      }
      return res;
    }
    ftype cov(vec const & a, vec const & b) {
      ftype ma = mean(a);
      ftype mb = mean(b);
      ftype res = 0.0;
      for( int ia = 0; ia < a.size(); ++ia) {
        for( int ib = 0; ib < b.size(); ++ib) {
          res += (ma-a(ia))*(mb-b(ib));
        }
      }
      return res;
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
