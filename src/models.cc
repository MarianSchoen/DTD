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
    mat estimate_c_direct(mat const & x, mat const & y, vec const & g, mat const & xtgxi) {
      return xtgxi*x.transpose()*g.asDiagonal()*y;
    }
    ftype GoertlerModel::eval(vec const & g) const {
      double res(0.0);
      mat xtgxi = invxtgx(m_x, g);
      mat c_hat = estimate_c_direct(m_x,m_y,g,xtgxi);

      return 0.0;
    }
    void GoertlerModel::grad(vec & gr, vec const & g) const {
    }
  }
}
