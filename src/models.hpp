#pragma once
#include "config.hpp"
#include <limits>

namespace dtd {
  namespace models {
    enum class NormFunctions {
                              IDENTITY,
                              NORM2,
                              NORM1
    };
    enum class ThresholdFunctions {
                                   SOFTMAX
    };
    enum class SubspaceFunctions {
                                   POSITIVE
    };
    enum class CoeffEstimation {
                                DIRECT,
                                NNLS
    };
    mat invxtgx(mat const & x, vec const & g);
    class GoertlerModel {
    private:
      int m_ngenes, m_ncells, m_nsamples;
      mat m_x, m_y, m_c, m_xtgxi;
      vec m_g;
      NormFunctions m_normfn;
      ThresholdFunctions m_threshfn;
      SubspaceFunctions m_subspfn;
      CoeffEstimation m_estim_c;
    public:
      GoertlerModel(mat x, mat y, mat c, vec g, NormFunctions fn = NormFunctions::IDENTITY, ThresholdFunctions thresh = ThresholdFunctions::SOFTMAX, SubspaceFunctions subsp = SubspaceFunctions::POSITIVE, CoeffEstimation estimC = CoeffEstimation::DIRECT) : m_ngenes(x.rows()), m_ncells(x.cols()), m_nsamples(y.cols()), m_x(x), m_y(y), m_c(c), m_g(g), m_normfn(fn), m_threshfn(thresh), m_subspfn(subsp), m_estim_c(estimC) {
        // some error checking:
        if( m_x.rows() != m_ngenes )
          throw std::runtime_error("number of rows of X is not ngenes.");
        if( m_x.cols() != m_ncells)
          throw std::runtime_error("number of cols of X is not ncells.");
        if( m_y.rows() != m_ngenes)
          throw std::runtime_error("number of rows of Y is not ngenes.");
        if( m_y.cols() != m_nsamples)
          throw std::runtime_error("number of cols of Y is not nsamples.");
        if( m_c.rows() != m_ncells)
          throw std::runtime_error("number of rows of C is not ncells.");
        if( m_c.cols() != m_nsamples)
          throw std::runtime_error("number of cols of C is not nsamples.");
        if( m_g.size() != m_ngenes && m_g.cols() == 1)
          throw std::runtime_error("g is not a ngenes long vector.");

        m_xtgxi = invxtgx(m_x, m_g);
      }
      mat estimate_c(vec const & g) const;
      mat estimate_c(vec const & g, mat const & xtgxi) const;
      inline vec const & getParams() const { return m_g; }
      inline void setParams(vec const & g) {
        m_g = g;
        m_xtgxi = invxtgx(m_x, m_g);
      }
      ftype evaluate(vec const & params) const {
        return evaluate(params, invxtgx(m_x, params));
      }
      ftype evaluate(vec const & params, mat const & xtgxi) const;
      inline ftype evaluate() const { return evaluate(m_g, m_xtgxi); }
      inline void grad(vec & gr) const { grad_explicit_inverse(gr, m_g, m_xtgxi); }
      void grad(vec & gr, vec const & g) const;
      void grad_explicit_inverse(vec & gr, vec const & param, mat const & xtgxi) const;
      inline std::size_t dim() const {
        return static_cast<std::size_t>(m_ngenes);
      }
      vec threshold(vec const & v, ftype softfactor) const;
      void norm_constraint(vec & v) const;
      vec subspace_constraint(vec const & v) const;
    };
  }
}
