#include <set>
#include "Eigen/Core"
#include "Eigen/LU"
#include "config.hpp"
#include <limits>

void set_s( vec & s, mat const & a, std::set<int> const & r, std::set<int> const & p, vec const & y){
  mat ap = mat(y.size(), p.size());
  vec sp = vec(p.size());
  int i = 0;
  for( auto e : p ){
    ap.col(i++) = a.col(e);
  }
  auto ata = ap.transpose()*ap;
  sp = ata.inverse()*ap.transpose()*y;
  i = 0;
  for( auto e : p ){
    s(e) = sp(i++);
  }
  for( auto e : r ){
    s(e) = 0.0;
  }
}
ftype minimum( vec const & v, std::set<int> q ) {
  ftype m = v(*q.begin());
  for( auto e : q ){
    if( m > v(e) ) m = v(e);
  }
  return m;
}
vec nnls(mat const & a, vec const & y, ftype eps=1e3*std::numeric_limits<ftype>::epsilon()) {
  const auto m = a.rows();
  const auto n = a.cols();
  if( y.size() != m ) {
    throw std::runtime_error("nnls: incompatible input lengths of a and y.");
  }
  std::set<int> p, r;
  for( int i = 0; i < n; ++i ){
    r.insert(i);
  }
  vec x = vec::Zero(n), s(n);
  vec w = a.transpose()*(y - a*x);

  Eigen::Index j;
  while( (! r.empty() ) && w.maxCoeff(&j) > eps ) {
    if( r.count(j) == 0 )
      throw std::runtime_error("nnls: maximum of w is not in R.");
    p.insert(j);
    r.erase(j);

    set_s(s, a, r, p, y);

    while( (! p.empty()) && (minimum(s, p) <= 0.0) ) {
      std::set<ftype> alphas;
      for( auto e : p ){
        if( s[e] <= eps ) {
          alphas.insert(x[e] / (x[e] - s[e]));
        }
      }
      ftype alpha = *std::min_element(alphas.cbegin(), alphas.cend());
      x = x + alpha * (s - x);
      for( auto e : p ){
        if( x[e] == 0.0 ){
          r.insert(e);
          p.erase(e);
        }
      }
      set_s(s, a, r, p, y);
    }
    x = s;
    w = a.transpose()*(y - a*x);
  }
  return x;
}

