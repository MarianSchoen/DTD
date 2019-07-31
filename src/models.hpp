#pragma once
#include "config.hpp"
#include <limits>

namespace dtd {
  namespace models {
    enum class NormFunctions {
                              IDENTITY,
                              NORM2
    };
    enum class ThresholdFunctions {
                                   SOFTMAX
    };
    enum class SubspaceFunctions {
                                   POSITIVE
    };
    mat invxtgx(mat const & x, vec const & g);
    bool check_posdefmat(mat const & x, vec const & g, mat const & xtgxi, ftype eps = std::numeric_limits<ftype>::epsilon());
    class GoertlerModel {
    private:
      int m_ngenes, m_ncells, m_nsamples;
      mat m_x, m_y, m_c, m_xtgxi;
      vec m_g;
      bool m_healthy;
      NormFunctions m_normfn;
      ThresholdFunctions m_threshfn;
      SubspaceFunctions m_subspfn;
    public:
      GoertlerModel(mat x, mat y, mat c, vec g, NormFunctions fn = NormFunctions::IDENTITY, ThresholdFunctions thresh = ThresholdFunctions::SOFTMAX, SubspaceFunctions subsp = SubspaceFunctions::POSITIVE) : m_ngenes(x.rows()), m_ncells(x.cols()), m_nsamples(y.cols()), m_x(x), m_y(y), m_c(c), m_g(g), m_healthy(true), m_normfn(fn), m_threshfn(thresh), m_subspfn(subsp) {
        assert(m_x.rows() == m_ngenes);
        assert(m_x.cols() == m_ncells);
        assert(m_y.rows() == m_ngenes);
        assert(m_y.cols() == m_nsamples);
        assert(m_c.rows() == m_ncells);
        assert(m_c.cols() == m_nsamples);
        assert(m_g.size() == m_ngenes && m_g.cols() == 1);
        m_xtgxi = invxtgx(m_x, m_g);
        m_healthy &= check_posdefmat(m_x, m_g, m_xtgxi);
      }
      inline vec const & getParams() const { return m_g; }
      inline void setParams(vec const & g) {
        m_g = g;
        m_xtgxi = invxtgx(m_x, m_g);
        m_healthy &= check_posdefmat(m_x, m_g, m_xtgxi);
      }
      bool isHealthy() const { return m_healthy; }
      ftype evaluate(vec const & params) const {
        // even if this inversion is unsafe, the model itself will still remain "healthy", meaning the combination of g and x lead to a proper inversion.
        const mat xtgxi = invxtgx(m_x, params);
        return evaluate(params, xtgxi);
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
