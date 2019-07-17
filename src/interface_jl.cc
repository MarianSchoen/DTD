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
