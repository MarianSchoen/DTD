#include "Eigen/Core"
#include <functional>
#include <vector>
#include <iostream> //TODO: remove

#include "config.hpp"
#include "utils.hpp"

using namespace Eigen;

namespace dtd {
  namespace solvers {
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

    template<class Model>
    double FistaSolver<Model>::iterate() {
      std::cout << "hello from FistaSolver::iterate()\n";
      return 0.0;
    }
  }
}
