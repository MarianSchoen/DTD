#include "Eigen/Core"
#include <functional>
#include <vector>

// TODO: think if single precision is enough? probably.. or first solve on floats, then iterate on double?
typedef double ftype;

// typedef Eigen::Array<ftype, Eigen::Dynamic, 1> vec;
// TODO: for the reference matrix, a sparse matrix type may also be interesting?? Does it pay off?
// typedef Eigen::Array<ftype, Eigen::Dynamic, 2> mat;
using vec = Eigen::VectorXd;
using mat = Eigen::MatrixXd;


using namespace Eigen;

namespace dtd {
  mat invxtgx(mat const & x, VectorXd const & g);
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
  namespace stat {
    ftype mean(vec const & v);
    ftype var(vec const & v);
    ftype cov(vec const & a, vec const & b);
  }
  // A model m must implement:
  // - a loss function, ftype m.eval(vec const & g) const
  // - a gradient, void m.grad(vec & gr, vec const & g) const
  // - a function std::size_t m.dim() const returning the dimension (i.e. the size of g and gr above)
  template<class Model>
  class FistaSolver {
  private:
    vec m_grad, m_g;
    Model m_model;
  public:
    FistaSolver(Model m) : m_model(m) {
      m_grad.resize(m.dim());
      m_g.resize(m.dim());
    }
    void setG(vec const & g) { m_g = g; }
    vec getG() const { return m_g; }
    void updateGradient() {
      m_model.grad(m_grad, m_g);
    }
    // for external use:
    ftype feval() const {
      return m_model.eval(m_g);
    }
    ftype iterate(); // returns value of loss function
  };
}

// helper for external bindings to julia / R / python / ...
// TODO: there are certainly better ways for a plain copy.
inline vec vecFromPtr(ftype* d, unsigned int len) {
  vec vec(len);
  for( unsigned int i = 0; i < len; ++i){
    vec(i) = d[i];
  }
  return vec;
}

inline mat matFromPtr(ftype* d, unsigned int rows, unsigned int cols) {
  mat m(rows, cols); // TODO: check row major
  for( unsigned int i = 0; i < cols; ++i){
    for( unsigned int j = 0; j < rows; ++j)
      m(j,i) = d[i*rows+j];
  }
  return m;
}
