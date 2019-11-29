#pragma once
#include "config.hpp"

namespace dtd {
  namespace stat {
    ftype var(vec const & v);
    ftype std(vec const & a);
    ftype cov(vec const & a, vec const & b);
    ftype cor(vec const & a, vec const & b);
  }
}
