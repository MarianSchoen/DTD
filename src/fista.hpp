#pragma once
#include "Eigen/Core"
#include <functional>
#include <vector>

#include "config.hpp"
#include "utils.hpp"

using namespace Eigen;

namespace dtd {
  namespace solvers {
    class AvgResidualChecker {
    private: 
      unsigned int m_navg, m_iteration;
      ftype m_maxres, m_sum;
      std::vector<ftype> m_deltas;
    public:
      AvgResidualChecker(unsigned int navg, ftype maximum_residual)
          : m_navg(navg), m_iteration(0), m_maxres(maximum_residual),
            m_sum(0.0), m_deltas(navg, 0.0) {}
      bool operator()(ftype delta_loss) {
        auto index = m_iteration%m_navg;
        m_sum -= m_deltas[index];
        m_deltas[index] = delta_loss;
        m_sum += delta_loss;
        m_iteration += 1;
        return ( m_iteration >= m_navg ) && ( m_sum / m_navg < m_maxres );
      }
    };
    // Barzilai & Borwein (1988) way of determining step sizes:
    template<class Model>
    double bb_learning_rate(Model const & m, vec const & params) {
      // requires two gradient evaluations
      vec grad, grad2;
      m.grad(grad, params);
      vec deltaX = - grad / grad.norm();
      m.grad(grad2, params + deltaX);
      vec deltaY = grad2 - grad;
      return abs(deltaX.dot(deltaY)) / deltaY.dot(deltaY);
    }

    inline double nesterov_factor(int k) {
      return 2.0 / static_cast<double>(k+1);
    }

    // A model m must implement:
    // - a loss function, ftype m.evaluate(vec const & g) const
    // - a gradient, void m.grad(vec & gr, vec const & g) const
    // - a function std::size_t m.dim() const returning the dimension (i.e. the size of g and gr above)
    template<class Model>
    class FistaSolver {
    private:
      vec m_grad;
      ftype m_learning_rate, m_linesearch_speed;
      bool m_restarts;
      int m_cyclelength, m_nesterov_counter;
      ftype m_delta_y_1, m_delta_y_2;
    public:
      FistaSolver(ftype learningrate, ftype linesearchspeed, int cycles, bool restarts) : m_nesterov_counter(2), m_delta_y_1(0.0), m_delta_y_2(0.0)
      {
        setLearningRate(learningrate);
        setLinesearchSpeed(linesearchspeed);
        setCyclelength(cycles);
        enableRestart(restarts);
      }
      // ctor with some defaults:
      FistaSolver() : FistaSolver(0.1, 2.0, 5, true) {}

      // algorithm parameters:
      void setLearningRate(ftype learningrate) {
        if( learningrate <= 0 )
          throw std::runtime_error("non-positive learning rate.");
        m_learning_rate = learningrate;
      }
      int getNesterovCounter() const {
        return m_nesterov_counter;
      }
      ftype getDeltaYBeforeNesterov() const {
        return m_delta_y_1;
      }
      ftype getDeltaYAfterNesterov() const {
        return m_delta_y_2;
      }
      ftype getLearningRate() const { return m_learning_rate; }
      void setLearningAuto(Model const & m) {
        setLearningRate(bb_learning_rate(m, m.getParams()));
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
      // for external use:
      ftype feval(Model const & m) const {
        return m.evaluate();
      }
      int solve(Model & m, int iter, double epsilon, std::size_t navg, double lambda, std::function<void(Model const & m, vec const &)> callback);
      inline int solve(Model & m, int iter, double epsilon, std::size_t navg, double lambda) { return solve(m, iter, epsilon, lambda, [](Model const &, vec const &){}); }
    };

    template<class Model>
    int FistaSolver<Model>::solve(Model & model, int maxiter, double epsilon, std::size_t navg, double lambda, std::function<void(Model const &, vec const &)> callback) {
      auto is_converged = AvgResidualChecker(navg, epsilon);
      vec y_vec = model.getParams();
      model.norm_constraint(y_vec);
      vec g_new = y_vec;
      ftype fy = model.evaluate();
      ftype fy_old = fy;

      mat testmat = mat(m_cyclelength, model.dim());
      vec testy   = vec(m_cyclelength);
      vec nesterov_dir = vec(model.dim());

      int iter = 2;
      for( ; iter <= maxiter; ++iter ){
        ftype rate = m_learning_rate / static_cast<double>(m_cyclelength - 1);
        model.grad(m_grad, y_vec);
        for(int i = 0; i < m_cyclelength; ++i) {
          ftype alpha = rate * static_cast<ftype>(i);
          testmat.row(i) = model.threshold(y_vec - alpha * m_grad, alpha * lambda);
          testy(i) = model.evaluate(testmat.row(i));
        }
        int minindex;
        ftype fy = testy.minCoeff(&minindex);

        if( minindex == 0 ) {
          // undershoot -> decrease step size
          m_learning_rate /= m_linesearch_speed;
        } else if ( minindex == m_cyclelength - 1 ) {
          // overshoot -> increase step size
          m_learning_rate *= m_linesearch_speed;
        }


        if( fy_old > fy ) {
          // descent
          y_vec = g_new;
          g_new = testmat.row(minindex);
          model.norm_constraint(g_new);
          m_nesterov_counter += 1;
        } else {
          // ascent
          // reset f_y to its previous value
          model.setParams(g_new);
          fy = model.evaluate(g_new);
          if( m_restarts )
            m_nesterov_counter = 2;
        }

        m_delta_y_1 = fy_old - fy;

        callback(model, g_new);
        fy_old = fy;

        // linesearch for nesterov extrapolation
        ftype factor = nesterov_factor(m_nesterov_counter);
        ftype deltafac = factor / static_cast<double>(m_cyclelength - 1);
        nesterov_dir = g_new - y_vec; // 0 up to normalization of g_new
        for( int i = 0; i < m_cyclelength; ++i) {
          ftype alpha = static_cast<ftype>(i)*deltafac;
          testmat.row(i) = model.subspace_constraint(g_new + alpha * nesterov_dir);
          testy(i) = model.evaluate(testmat.row(i));
        }
        fy = testy.minCoeff(&minindex);
        m_delta_y_2 = fy_old - fy;
        y_vec = testmat.row(minindex);
        model.norm_constraint(y_vec);

        if( is_converged(m_delta_y_1 + m_delta_y_2) )
          break; // converged.
      }
      model.setParams(g_new);
      return iter;
    }
  }
}
