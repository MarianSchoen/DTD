#include "models.hpp"
#include "utils.hpp"
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
      assert(g.size() == m_ngenes);
      double res(0.0);
      mat xtgxi = invxtgx(m_x, g);
      assert(xtgxi.rows() == m_ncells);
      assert(xtgxi.cols() == m_ncells);
      mat c_hat = estimate_c_direct(m_x,m_y,g,xtgxi);
      assert(c_hat.rows() == m_ncells);
      assert(c_hat.cols() == m_nsamples);

      for( std::size_t icell = 0; icell < m_ncells; ++icell ){
        res -= stat::cor(m_c.row(icell), c_hat.row(icell));
      }
      return res;
    }
    void clampPos(vec & x) {
      for( unsigned int i = 0; i < x.size(); ++i){
        x.coeffRef(i) = std::min(x.coeff(i), 0.0);
      }
    }
    void GoertlerModel::grad(vec & gr, vec const & g) const {
      assert(g.size() == m_ngenes);
      if( gr.size() != g.size() )
        gr.resize(g.size());
      mat xtgxi = invxtgx(m_x, g);
      mat c_hat = estimate_c_direct(m_x,m_y,g,xtgxi);
      assert(c_hat.rows() == m_ncells);
      assert(c_hat.cols() == m_nsamples);

      // naively:
      mat A = mat(m_ncells, m_nsamples);
      for( unsigned int icell = 0; icell < m_ncells; ++icell){
        ftype std_c_hat  = stat::std(c_hat.row(icell));
        ftype std_c      = stat::std(m_c.row(icell));
        ftype mean_c_hat = c_hat.row(icell).mean();
        ftype mean_c     = m_c.row(icell).mean();
        A.row(icell) = (stat::cov(m_c.row(icell), c_hat.row(icell)) / (m_nsamples * std_c_hat * std_c_hat) * (c_hat.row(icell).array() - mean_c_hat) - (m_c.row(icell).array() - mean_c)/m_nsamples) / (std_c * std_c_hat);
      }
      // TODO: speed this up by ONLY computing the diagonal
      mat tmp = (m_y - m_x*xtgxi*m_x.transpose()*g.asDiagonal()*m_y)*A.transpose()*xtgxi*m_x.transpose();
      gr = tmp.diagonal();
      clampPos(gr);
    }
  }
}
