#include "utils.hpp"
#include <cmath>
#include <limits>
namespace dtd {
  namespace stat {
    ftype var(vec const & v) {
      vec vm = v.array() - v.mean();
      auto n = v.size();
      if( n <= 1 )
        throw std::runtime_error("cannot compute variance for sample sizes < 2.");
      return vm.dot(vm) / (n-1);
    }
    ftype std(vec const & a) {
      return std::sqrt(var(a));
    }
    ftype cov(vec const & a, vec const & b) {
      int n = a.size();
      assert(n == b.size());
      vec am = a.array() - a.mean();
      vec bm = b.array() - b.mean();
      return am.dot(bm)/(n-1);
    }
    ftype cor(vec const & a, vec const & b) {
      auto va = var(a);
      auto vb = var(b);
      if( va < a.size()*std::numeric_limits<ftype>::epsilon() ||
          vb < b.size()*std::numeric_limits<ftype>::epsilon() )
        throw std::runtime_error("cor: cannot compute correlation of things that don't vary.");
      assert(va > 0.0 && vb > 0.0);
      return cov(a, b) / std::sqrt(var(a)*var(b));
    }
  }
}
