#pragma once
#include "config.hpp"

namespace dtd {
  namespace models {
  mat invxtgx(mat const & x, vec const & g);
  class GoertlerModel {
  private:
    int m_ngenes;
    mat m_refmat, m_measdat, m_c;
  public:
    ftype eval(vec const & g) const;
    void grad(vec & gr, vec const & g) const;
    std::size_t dim() const {
      return static_cast<std::size_t>(m_ngenes);
    }
  };
}
}
