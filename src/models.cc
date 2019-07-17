#include "models.hpp"
#include <Eigen/Cholesky>
#include <iostream> // TODO: remove

namespace dtd {
  namespace models {
    mat invxtgx(mat const & x, vec const & g) {
      std::cout << "x: " << x.rows() << " x " << x.cols() << " matrix.\n";
      std::cout << "g: " << g.size() << " vector.\n";
      mat xtgx = x.transpose()*g.asDiagonal()*x;
      mat xi = xtgx.llt().solve(mat::Identity(x.cols(), x.cols()));
      std::cout <<"from c++:\n"<< xi << "\n";
      return xi;
    }
    ftype GoertlerModel::eval(vec const & g) const {
      return 0.0;
    }
    void GoertlerModel::grad(vec & gr, vec const & g) const {
    }
  }
}
