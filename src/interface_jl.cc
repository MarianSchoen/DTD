#include "dtd.hpp"
#include "Eigen/Core"
#include "interface_jl.hpp"
Eigen::MatrixXd vecToEMat(std::vector<double> const & mat, unsigned int rows, unsigned int cols) {
  if( cols*rows != mat.size() )
    throw std::runtime_error("vector and matrix dimensions don't match.");

  MatrixXd res( rows, cols);

  for( unsigned int i = 0; i < cols; ++i ){
    for( unsigned int j = 0; j < rows; ++j ){
      res(j, i) = mat[i*rows+j];
    }
  }
  return res;
}
Eigen::VectorXd vecToEVec(std::vector<double> const & vec) {
  VectorXd res(vec.size());
  for( unsigned int i = 0; i < vec.size(); ++i) {
    res(i) = vec[i];
  }
  return res;
}
std::vector<double> ematToVec(Eigen::MatrixXd const & m) {
  std::vector<double> res(m.size());
  for( unsigned int i = 0; i < m.cols(); ++i) {
    for( unsigned int j = 0; j < m.rows(); ++j) {
      res[i*m.rows()+j] = m(j,i);
    }
  }
  return res;
}
std::vector<double> evecToVec(Eigen::VectorXd const & v) {
  std::vector<double> res(v.size());
  for( unsigned int i = 0; i < v.size(); ++i)
    res[i] = v[i];
  return res;
}

std::vector<double> ptrToVec(double const * ptr, unsigned int len) {
  std::vector<double> res(len);
  for( unsigned int i = 0; i < len; ++i){
    res[i] = ptr[i];
  }
  return res;
}
std::vector<double> xtgxinv(std::vector<double> const & xv, std::vector<double> const & gv, std::size_t ngenes, std::size_t ncells) {
  MatrixXd x = vecToEMat(xv, ngenes, ncells);
  MatrixXd g = vecToEVec(gv);
  return ematToVec(dtd::models::invxtgx(x,g));
}

double var(std::vector<double> const & x) {
  return dtd::stat::var(vecToEVec(x));
}
double cov(std::vector<double> const & x, std::vector<double> const & y) {
  return dtd::stat::cov(vecToEVec(x), vecToEVec(y));
}

dtd::models::GoertlerModel makeGoertlerModel(std::vector<double> const & xv, std::vector<double> const & yv, std::vector<double> const & cv, std::vector<double> const & gv, std::size_t ngenes, std::size_t ncells, std::size_t nsamples) {
  MatrixXd x = vecToEMat(xv, ngenes, ncells);
  MatrixXd y = vecToEMat(yv, ngenes, nsamples);
  MatrixXd c = vecToEMat(cv, ncells, nsamples);
  VectorXd g = vecToEVec(gv);
  assert(g.size() == ngenes);
  return dtd::models::GoertlerModel(x,y,c,g);
}

double evalGoertlerModel(dtd::models::GoertlerModel const * model, std::vector<double> const & g) {
  return model->evaluate(vecToEVec(g));
}

std::vector<double> gradGoertlerModel(dtd::models::GoertlerModel const * model, std::vector<double> const & g) {
  VectorXd grad(g.size());
  model->grad(grad, vecToEVec(g));
  return ematToVec(grad);
}

double solveFista(dtd::models::GoertlerModel * model, std::vector<double> & params, double lambda, std::size_t maxiter){
  dtd::solvers::FistaSolver<dtd::models::GoertlerModel> solver; // TODO: pass algorithmic params...
  solver.solve(*model, maxiter, lambda);
  params = evecToVec(model->getParams());
  return model->evaluate();
}

double bb_learning_rate(dtd::models::GoertlerModel const * model, std::vector<double> const & params) {
  // TODO:  interface the other parameters of the fista solver, e.g., lambda, linesearchspeed, etc. pp.
  // what is a good way to do this?
  vec g = vecToEVec(params);
  return dtd::solvers::bb_learning_rate(*model, g);
}
