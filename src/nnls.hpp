#include <set>
#include "Eigen/Core"
#include "config.hpp"
#include <limits>

void set_s( vec & s, mat const & a, std::set<int> const & r, std::set<int> const & p, vec const & y);
ftype minimum( vec const & v, std::set<int> q );
vec nnls(mat const & a, vec const & y, ftype eps=1e3*std::numeric_limits<ftype>::epsilon());
