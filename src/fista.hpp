#include "Eigen/Core"
#include <functional>
#include <vector>

#include "config.hpp"
#include "utils.hpp"

using namespace Eigen;

namespace dtd {
  namespace solvers {
    // Barzilai & Borwein (1988) way of determining step sizes:
    template<class Model>
    double bb_learning_rate(Model const & m, vec const & params) {
      // requires two gradient evaluations
      vec grad, grad2;
      m.grad(grad, params);
      vec deltaX = - grad / grad.norm();
      m.grad(grad2, params + deltaX);
      vec deltaY = grad2 - grad;
      return deltaX.dot(deltaY) / deltaY.norm();
    }

    inline double nesterov_factor(int k) {
      return 2.0 / static_cast<double>(k+1);
    }

    // A model m must implement:
    // - a loss function, ftype m.eval(vec const & g) const
    // - a gradient, void m.grad(vec & gr, vec const & g) const
    // - a function std::size_t m.dim() const returning the dimension (i.e. the size of g and gr above)
    template<class Model>
    class FistaSolver {
    private:
      vec m_grad, m_g;
      ftype m_learning_rate, m_linesearch_speed;
      bool m_restarts;
      int m_cyclelength, m_nesterov_counter;
      Model m_model;
    public:
      FistaSolver(Model m, vec g, ftype learningrate, ftype linesearchspeed, int cycles, bool restarts) : m_model(m), m_g(g), m_nesterov_counter(2)
      { // TODO: avoid copy of model? keep it as a const ref?
        assert(m.dim() == m_g.size());
        m_grad.resize(m.dim());
        setLearningRate(learningrate);
        setLinesearchSpeed(linesearchspeed);
        setCyclelength(cycles);
        enableRestart(restarts);
      }
      FistaSolver(Model m, ftype learningrate, ftype linesearchspeed, int cycles, bool restarts) : FistaSolver(m, vec(m.dim()), learningrate, linesearchspeed, cycles, restarts)
      {}
      // ctor with some defaults:
      FistaSolver(Model m) : FistaSolver(m, 0.1, 2.0, 5, true) {}

      // algorithm parameters:
      void setLearningRate(ftype learningrate) {
        if( learningrate <= 0 )
          throw std::runtime_error("non-positive learning rate.");
        m_learning_rate = learningrate;
      }
      ftype getLearningRate() const { return m_learning_rate; }
      void setLearningAuto() {
        setLearningRate(bb_learning_rate(m_model, m_g));
      }
      void setLinesearchSpeed(ftype lsspeed) {
        if( lsspeed <= 1 )
          throw std::runtime_error("linesearch speed must be > 1.");
        m_linesearch_speed = lsspeed;
      }
      ftype getLinesearchSpeed() const { return m_linesearch_speed; }
      void enableRestart(bool rs = true) { m_restarts = rs; }
      ftype isRestarting() const { return m_restarts; }
      void setCyclelength( int cyclelen) {
        if( cyclelen <= 0 )
          throw std::runtime_error("non-positive cycle length!");
        m_cyclelength = cyclelen;
      }
      int getCyclelength() const { return m_cyclelength; }
      void setG(vec const g) { m_g = g; }
      vec const & getG() const { return m_g; }
      // for external use:
      ftype feval() const {
        return m_model.eval(m_g);
      }
      void solve(std::size_t iter, double lambda); // returns value of loss function
    };

    template<class Model>
    void FistaSolver<Model>::solve(std::size_t maxiter, double lambda) {
      vec y_vec = m_g;
      vec g_new = m_g;
      vec u_vec, g_old;
      ftype fy = m_model.eval(m_g);
      ftype fy_old = fy;

      // TODO preallocate dynamic memory

      for( std::size_t iter = 2; iter <= maxiter; ++iter ){
        ftype rate = m_learning_rate / static_cast<double>(m_cyclelength - 1);
        m_model.grad(m_grad, y_vec);
        // TODO if memory is an issue, rethink this:
        mat testmat = mat(m_cyclelength, m_model.dim());
        vec testy   = vec(m_cyclelength);
        for(int i = 0; i < m_cyclelength; ++i) {
          double alpha = rate * i;
          testmat.row(i) = m_model.threshold(y_vec - alpha * m_grad, alpha * lambda);
          testy(i) = m_model.eval(testmat.row(i));
        }
        std::size_t minindex;
        ftype fy = testy.minCoeff(&minindex);

        if( minindex == 0 ) {
          // undershoot -> decrease step size
          m_learning_rate /= m_linesearch_speed;
        } else if ( minindex == m_cyclelength - 1 ) {
          // overshoot -> increase step size
          m_learning_rate *= m_linesearch_speed;
        }

        u_vec = testmat.row(minindex);
        m_model.norm_constraint(u_vec);

        g_old = g_new;

        if( fy_old > fy ) {
          // descent
          // TODO: swap instead?! (here and below)
          g_new = u_vec;
        } else {
          // ascent
          fy = m_model.eval(g_new);
          if( m_restarts )
            m_nesterov_counter = 2;
        }

        fy_old = fy;

        // linesearch for nesterov extrapolation
        ftype factor = nesterov_factor(m_nesterov_counter);
        ftype deltafac = factor / static_cast<double>(m_cyclelength - 1);
        vec nesterov_dir = g_new - g_old;
        for( int i = 0; i < m_cyclelength; ++i) {
          ftype alpha = i*deltafac;
          testmat.row(i) = m_model.subspace_constraint(g_new + alpha * nesterov_dir);
          testy(i) = m_model.eval(testmat.row(i));
        }
        testy.minCoeff(&minindex);
        y_vec = testmat.row(minindex);
      }
      m_model.norm_constraint(g_new);
      m_g.swap(g_new);
    }
  }
}
