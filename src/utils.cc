#include "utils.hpp"
namespace dtd {
  namespace stat {
    ftype var(vec const & v) {
      vec vm = v.array() - v.mean();
      auto n = v.size();
      return vm.dot(vm) / (n-1);
    }
    ftype cov(vec const & a, vec const & b) {
      int n = a.size();
      assert(n == b.size());
      vec am = a.array() - a.mean();
      vec bm = b.array() - b.mean();
      return am.dot(bm)/(n-1);
    }
  }
}
