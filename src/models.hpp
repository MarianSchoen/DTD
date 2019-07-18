#pragma once
#include "config.hpp"

namespace dtd {
  namespace models {
    mat invxtgx(mat const & x, vec const & g);
    class GoertlerModel {
    private:
      int m_ngenes, m_ncells, m_nsamples;
      mat m_x, m_y, m_c;
    public:
      GoertlerModel(mat x, mat y, mat c) : m_x(x), m_y(y), m_c(c), m_ngenes(x.rows()), m_ncells(x.cols()), m_nsamples(y.cols()) {
        assert(m_x.rows() == m_ngenes);
        assert(m_x.cols() == m_ncells);
        assert(m_y.rows() == m_ngenes);
        assert(m_y.cols() == m_nsamples);
        assert(m_c.rows() == m_ncells);
        assert(m_c.cols() == m_nsamples);
      }
      ftype eval(vec const & g) const;
      void grad(vec & gr, vec const & g) const;
      std::size_t dim() const {
        return static_cast<std::size_t>(m_ngenes);
      }
    };
  }
}
