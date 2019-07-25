#pragma once
#include "config.hpp"

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
    class GoertlerModel {
    private:
      int m_ngenes, m_ncells, m_nsamples;
      mat m_x, m_y, m_c;
      vec m_g;
      NormFunctions m_normfn;
      ThresholdFunctions m_threshfn;
      SubspaceFunctions m_subspfn;
    public:
      GoertlerModel(mat x, mat y, mat c, vec g, NormFunctions fn = NormFunctions::IDENTITY, ThresholdFunctions thresh = ThresholdFunctions::SOFTMAX, SubspaceFunctions subsp = SubspaceFunctions::POSITIVE) : m_ngenes(x.rows()), m_ncells(x.cols()), m_nsamples(y.cols()), m_x(x), m_y(y), m_c(c), m_g(g), m_normfn(fn), m_threshfn(thresh), m_subspfn(subsp) {
        assert(m_x.rows() == m_ngenes);
        assert(m_x.cols() == m_ncells);
        assert(m_y.rows() == m_ngenes);
        assert(m_y.cols() == m_nsamples);
        assert(m_c.rows() == m_ncells);
        assert(m_c.cols() == m_nsamples);
        assert(m_g.size() == m_ngenes && m_g.cols() == 1);
      }
      inline vec const & getParams() const { return m_g; }
      void setParams(vec const & g) { m_g = g; }
      ftype evaluate(vec const & params) const;
      inline ftype evaluate() const { return evaluate(m_g); }
      void grad(vec & gr, vec const & g) const;
      void grad(vec & gr) const;
      inline std::size_t dim() const {
        return static_cast<std::size_t>(m_ngenes);
      }
      vec threshold(vec const & v, ftype softfactor) const;
      void norm_constraint(vec & v) const;
      vec subspace_constraint(vec const & v) const;
    };
  }
}
