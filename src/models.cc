#include "models.hpp"
#include "utils.hpp"
#include <Eigen/Cholesky>
#include <limits>
#include <iostream>

namespace dtd {
  namespace models {
    template<typename T>
    int sgn(T x) {
      return (T(0) < x) - (x < T(0));
    }
    vec softmax(vec const & v, ftype softfactor) {
      vec res(v.size());
      // TODO: non-loop based approach may be faster.
      for( int i = 0; i < v.size(); ++i ){
        double const & x = v.coeff(i);
        res.coeffRef(i) = sgn(x)*std::max(std::abs(x)-softfactor, 0.0);
      }
      return res;
    }
    mat invxtgx(mat const & x, vec const & g) {
      if( g.minCoeff() < 0 )
        throw std::runtime_error("invxtgx: g has negative entries and thus, x^T G x is not positive semi-definite and cannot be inverted.");
      mat xtgx = x.transpose()*g.asDiagonal()*x;
      // inversion happens here:
      mat xi = xtgx.llt().solve(mat::Identity(x.cols(), x.cols()));
      // the following is just a check. it is a bit costly, but not overly.
      const std::size_t n = x.cols();
      const ftype eps = std::numeric_limits<ftype>::epsilon();
      const mat zero = xi*xtgx - mat::Identity(n, n);
      if( zero.norm() > n*n*eps )
        throw std::runtime_error("invxtgx: could not invert x^T G x, probably because the rank of x is too low to compensate for too many zeros in g.");
      return xi;
    }
    mat estimate_c_direct(mat const & x, mat const & y, vec const & g, mat const & xtgxi) {
      return xtgxi*x.transpose()*g.asDiagonal()*y;
    }
    ftype GoertlerModel::evaluate(vec const & g, mat const & xtgxi) const {
      double res(0.0);
      assert(g.size() == m_ngenes);
      assert(xtgxi.rows() == m_ncells);
      assert(xtgxi.cols() == m_ncells);
      mat c_hat = estimate_c_direct(m_x,m_y,g,xtgxi);
      assert(c_hat.rows() == m_ncells);
      assert(c_hat.cols() == m_nsamples);

      try {
        for( int icell = 0; icell < m_ncells; ++icell ){
          res -= stat::cor(m_c.row(icell), c_hat.row(icell));
        }
      } catch(...) {
        throw;
      }
      return res / m_ncells;
    }
    inline void clampPos(vec & x) {
      for( unsigned int i = 0; i < x.size(); ++i){
        x.coeffRef(i) = std::min(x.coeff(i), 0.0);
      }
    }
    inline void clampNeg(vec & x) {
      for( unsigned int i = 0; i < x.size(); ++i){
        x.coeffRef(i) = std::max(x.coeff(i), 0.0);
      }
    }
    void GoertlerModel::grad_explicit_inverse(vec & gr, vec const & g, mat const & xtgxi) const {
      assert( g.size() == m_ngenes);
      if( gr.size() != m_ngenes )
        gr.resize(m_ngenes);
      mat c_hat = estimate_c_direct(m_x,m_y,g,xtgxi);
      assert(c_hat.rows() == m_ncells);
      assert(c_hat.cols() == m_nsamples);

      // naively:
      mat A = mat(m_ncells, m_nsamples);
      for( int icell = 0; icell < m_ncells; ++icell){
        ftype std_c_hat  = stat::std(c_hat.row(icell));
        ftype std_c      = stat::std(m_c.row(icell));
        ftype mean_c_hat = c_hat.row(icell).mean();
        ftype mean_c     = m_c.row(icell).mean();
        if( std_c_hat <= c_hat.rows()*std::numeric_limits<ftype>::epsilon() ||
            std_c     <= m_c.rows()*std::numeric_limits<ftype>::epsilon() )
          throw std::runtime_error("grad_explicit_inverse: cannot compute gradient, because there is no variance in C over the samples.");
        A.row(icell) = (stat::cov(m_c.row(icell), c_hat.row(icell)) / (m_nsamples * std_c_hat * std_c_hat) * (c_hat.row(icell).array() - mean_c_hat) - (m_c.row(icell).array() - mean_c)/m_nsamples) / (std_c * std_c_hat);
      }
      // TODO: speed this up by ONLY computing the diagonal
      mat tmp = (m_y - m_x*xtgxi*m_x.transpose()*g.asDiagonal()*m_y)*A.transpose()*xtgxi*m_x.transpose();
      gr = tmp.diagonal();
      clampPos(gr);
    }
    void GoertlerModel::grad(vec & gr, vec const & param) const {
      const mat xtgxi = invxtgx(m_x, param);
      grad_explicit_inverse(gr, param, xtgxi);
    }
    vec GoertlerModel::threshold(vec const & v, ftype softfactor) const {
      if( m_threshfn == ThresholdFunctions::SOFTMAX )
        return softmax(v, softfactor);
      else
        throw std::runtime_error("unimplemented threshold function.");
    }
    void GoertlerModel::norm_constraint(vec & v) const {
      if( m_normfn == NormFunctions::IDENTITY )
        ; // leave v unchanged
      else if( m_normfn == NormFunctions::NORM2 )
        v /= v.norm();
      else
        throw std::runtime_error("unimplemented norm function. ");
    }
    vec GoertlerModel::subspace_constraint(vec const & v) const {
      if( m_subspfn == SubspaceFunctions::POSITIVE ) {
        // set all negative values to zero:
        vec res(v);
        clampNeg(res);
        return res;
      }
      else
        throw std::runtime_error("unimplemented subspace constraint.");

    }
  }
}
