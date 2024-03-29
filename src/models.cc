#include "models.hpp"
#include "utils.hpp"
#include "nnls.hpp"
#include <Eigen/Cholesky>
#include <limits>
#include <sstream>

namespace dtd {
  namespace models {
    template<typename T>
    int sgn(T x) {
      return (T(0) < x) - (x < T(0));
    }
    vec softmax(vec const & v, ftype softfactor) {
      vec res(v.size());
      for( int i = 0; i < v.size(); ++i ){
        double const & x = v.coeff(i);
        res.coeffRef(i) = sgn(x)*std::max(std::abs(x)-softfactor, 0.0);
      }
      return res;
    }
    mat invxtgx(mat const & x, vec const & g, ftype eps = 1000*std::numeric_limits<ftype>::epsilon()) {
      if( g.minCoeff() < 0 )
        throw std::runtime_error("invxtgx: g has negative entries and thus, x^T G x is not positive semi-definite and cannot be inverted.");
      mat xtgx = x.transpose()*g.asDiagonal()*x;
      // inversion happens here, using Cholesky decomposition of the positive definite matrix xtgx:
      mat xi = xtgx.llt().solve(mat::Identity(x.cols(), x.cols()));
      // the following is just a check. it is a bit costly, but not overly.
      const std::size_t n = x.cols();
      // there is some empirical arbitrariness in here...
      // but this is probably still save (whenever it failed it was O(1) and not O(1e-13)
      const mat zero = xi*xtgx - mat::Identity(n, n);
      if( zero.trace() > n*eps ){
        std::stringstream ss;
        ss << "invxtgx: could not invert X^T diag(g) X.\n";
        ss << "Usually this happens when the rank of X is too low to compensate for too many zeros in g\n";
        ss << "It may help to: \n";
        ss << " * decrease lambda \n";
        ss << " * increase the number of features \n";
        ss << " * increase number of lambdas\n\n";
        ss << "in 'train_deconvolution_model' set 'cv.verbose=TRUE' for logging information\n";
        ss << "residual: " << zero.trace() << "\n";
        ss << "admissible threshold: " << eps << "\n";
        throw std::runtime_error(ss.str());
      }
      return xi;
    }
    mat estimate_c_direct(mat const & x, mat const & y, vec const & g, mat const & xtgxi) {
      return xtgxi*x.transpose()*g.asDiagonal()*y;
    }
    mat estimate_c_nnls(mat const & x, mat const & y, vec const & g, ftype eps = 1000*std::numeric_limits<ftype>::epsilon(), int maxiter = 10000) {
      mat gx = g.asDiagonal()*x;
      mat res = mat(x.cols(), y.cols());
      for( int i = 0; i < y.cols(); ++i ){
        res.col(i) = nnls(gx, g.asDiagonal()*y.col(i), eps, maxiter);
      }
      return res;
    }
    mat GoertlerModel::estimate_c(vec const & g) const {
      if( m_estim_c == CoeffEstimation::DIRECT ) {
        mat xtgxi = invxtgx(m_x, g, inv_prec);
        return estimate_c_direct(m_x, m_y, g, xtgxi);
      } else if( m_estim_c == CoeffEstimation::NNLS ){
        return estimate_c_nnls(m_x, m_y, g, inv_prec, 10000); // <- a large random number (maxiter)
      } else {
        throw std::runtime_error("unimplemented estimate_C function.");
      }
    }
    mat GoertlerModel::estimate_c(vec const & g, mat const & xtgxi) const {
      if( m_estim_c == CoeffEstimation::DIRECT ) {
        return estimate_c_direct(m_x, m_y, g, xtgxi);
      } else if( m_estim_c == CoeffEstimation::NNLS ){
        return estimate_c_nnls(m_x, m_y, g);
      } else {
        throw std::runtime_error("unimplemented estimate_C function.");
      }
    }
    ftype GoertlerModel::evaluate(vec const & g, mat const & xtgxi) const {
      double res(0.0);
      assert(g.size() == m_ngenes);
      assert(xtgxi.rows() == m_ncells);
      assert(xtgxi.cols() == m_ncells);
      mat c_hat = this->estimate_c(g);
      assert(c_hat.rows() == m_ncells);
      assert(c_hat.cols() == m_nsamples);

      try {
        for( int icell = 0; icell < m_ncells; ++icell ){
          if ( stat::std(c_hat.row(icell)) > c_hat.rows()*std::numeric_limits<ftype>::epsilon() )
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
      mat c_hat = this->estimate_c(g, xtgxi);
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
          // no variance in C, set sigma = 1, cov = 1:
          A.row(icell) = (c_hat.row(icell).array() - mean_c_hat) - (m_c.row(icell).array() - mean_c)/m_nsamples;
        else
          A.row(icell) = (stat::cov(m_c.row(icell), c_hat.row(icell)) / (m_nsamples * std_c_hat * std_c_hat) * (c_hat.row(icell).array() - mean_c_hat) - (m_c.row(icell).array() - mean_c)/m_nsamples) / (std_c * std_c_hat);
      }
      // TODO: speed this up by ONLY computing the diagonal
      mat tmp = (m_y - m_x*xtgxi*m_x.transpose()*g.asDiagonal()*m_y)*A.transpose()*xtgxi*m_x.transpose();
      gr = tmp.diagonal();
      clampPos(gr);
    }
    void GoertlerModel::grad(vec & gr, vec const & param) const {
      const mat xtgxi = invxtgx(m_x, param, inv_prec);
      grad_explicit_inverse(gr, param, xtgxi);
    }
    vec GoertlerModel::threshold(vec const & v, ftype softfactor) const {
      if( m_threshfn == ThresholdFunctions::SOFTMAX )
        return softmax(v, softfactor);
      else
        throw std::runtime_error("unimplemented threshold function.");
    }
    void GoertlerModel::norm_constraint(vec & v) const {
      ftype scale = 1.0;
      if( m_normfn == NormFunctions::IDENTITY ) {
        return; // leave v unchanged
      } else if( m_normfn == NormFunctions::NORM2 ) {
        //scale = v.size() / v.norm();
        scale = v.norm();
      } else if( m_normfn == NormFunctions::NORM1 ) {
        scale = v.lpNorm<1>();
      } else
        throw std::runtime_error("unimplemented norm function. ");
      v *= v.size() / scale;
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
