#include "dtd.hpp"
#include "Eigen/Core"
#include <R.h>
#include <Rdefines.h>

// copied from https://github.com/rehbergT/dgemm/blob/master/dgemmR/src/wrapper.cpp
SEXP getElementFromRList(SEXP RList, const char* name) {
    SEXP element = R_NilValue;
    SEXP names = getAttrib(RList, R_NamesSymbol);
    for (uint32_t i = 0; i < (uint32_t)length(RList); i++) {
        if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
            element = VECTOR_ELT(RList, i);
            break;
        }
    }
    return element;
}
MatrixXd getMatrixFromR(SEXP mat) {
  double* dmat = REAL(mat);
  int m = INTEGER(GET_DIM(mat))[0];
  int n = INTEGER(GET_DIM(mat))[1];
  MatrixXd res(m,n);
  for( int i = 0; i < m; ++i) {
    for( int j = 0; j < n; ++j) {
      res(i,j) = dmat[i*m+j];
    }
  }
  return res;
}

void fillPtr(double* dest, MatrixXd const & src) {
  // no bounds checking!!
  for( int i = 0; i < src.rows(); ++i ){
    for( int j = 0; j < src.cols(); ++j) {
      dest[i*src.rows()+j] = src(i,j);
    }
  }
}

SEXP dtd_solve_fista_goertler(SEXP tweak_, SEXP model_, SEXP _lambda, SEXP _maxiter){
  auto g = getMatrixFromR(tweak_);
  auto x = getMatrixFromR(getElementFromRList(model_, "X"));
  auto y = getMatrixFromR(getElementFromRList(model_, "Y"));
  auto c = getMatrixFromR(getElementFromRList(model_, "c"));
  double lambda = *REAL(_lambda);
  int maxiter = *INTEGER(_maxiter);

  dtd::models::GoertlerModel model(x,y,c);
  dtd::solvers::FistaSolver<dtd::models::GoertlerModel> solver(model);
  solver.setG(g);
  solver.solve( maxiter, lambda);

  SEXP result;
  PROTECT(result = allocMatrix(REALSXP, g.size(), 1));
  fillPtr(REAL(result), g);
  UNPROTECT(1);
  return result;
}
extern "C" {
  static const R_CallMethodDef callMethods[] = {
                                                { "dtd_solve_fista_goertler", (DL_FUNC)&dtd_solve_fista_goertler, 3},
                                                {NULL, NULL, 0}
  };
  void R_init_DTD(DllInfo* info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, TRUE);
  }
}
