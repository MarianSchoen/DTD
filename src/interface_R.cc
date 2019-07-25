#include "dtd.hpp"
#include "Eigen/Core"
#include <R.h>
#include <Rdefines.h>
#include <array>

// copied from https://github.com/rehbergT/dgemm/blob/master/dgemmR/src/wrapper.cpp
SEXP getElementFromRList(SEXP list_r, const char* name) {
    SEXP names = getAttrib(list_r, R_NamesSymbol);
    for (uint32_t i = 0; i < (uint32_t)length(list_r); i++) {
      Rprintf("DEBUG: cmp %s with %s...\n", CHAR(STRING_ELT(names, i)), name);
        if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
            Rprintf("DEBUG: element found..");
            return VECTOR_ELT(list_r, i);
        }
    }
    Rprintf("entry %s not found in list!\n", name);
    return R_NilValue;
}
MatrixXd getMatrixFromR(SEXP mat) {
  Rprintf("debug4a...\n");
  double* dmat = REAL(mat);
  Rprintf("debug4b...\n");
  int m = INTEGER(GET_DIM(mat))[0]; // TODO, see below !?!?!
  Rprintf("debug4c...\n");
  int n = INTEGER(GET_DIM(mat))[1];
  MatrixXd res(m,n);
  for( int i = 0; i < m; ++i) {
    for( int j = 0; j < n; ++j) {
      res(i,j) = dmat[i*m+j];
    }
  }
  return res;
}
VectorXd getVectorFromR(SEXP v) {
  Rprintf("debugv4a...\n");
  double* dv = REAL(v);
  Rprintf("debugv4b...\n");
  int m = XLENGTH(v);
  Rprintf("debugv4c...\n");
  vec res(m);
  for( int i = 0; i < m; ++i) {
      res(i) = dv[i];
  }
  return res;
}

void fillPtr(double* dest, MatrixXd const & src) {
  // no bounds checking!!
  for( int i = 0; i < src.rows(); ++i ){
    for( int j = 0; j < src.cols(); ++j) {
      dest[i*src.cols()+j] = src(i,j);
    }
  }
}

dtd::models::GoertlerModel make_model(SEXP model_) {
  auto x = getMatrixFromR(getElementFromRList(model_, "X"));
  auto y = getMatrixFromR(getElementFromRList(model_, "Y"));
  auto c = getMatrixFromR(getElementFromRList(model_, "C"));
  auto g = getVectorFromR(getElementFromRList(model_, "tweak"));

  return dtd::models::GoertlerModel(x,y,c, g);
}

SEXP dtd_solve_fista_goertler(SEXP model_, SEXP _lambda, SEXP _maxiter){
  double lambda = REAL(_lambda)[0];
  int maxiter = REAL(_maxiter)[0]; // TODO: somehow, integers are doubles, actually??
  auto model = make_model(model_);

  dtd::solvers::FistaSolver<dtd::models::GoertlerModel> solver(model);

  auto conv_vec = solver.solve( maxiter, lambda);

  const std::array<const char* const, 3> listnames = {"Tweak", "Convergence", "Lambda"};
  const std::size_t listlen = listnames.size();
  Rprintf("DEBUG: 10\n");
  SEXP listentrynames = PROTECT(allocVector(VECSXP, listlen));
  Rprintf("DEBUG: 11\n");
  for( auto i = 0u; i < listlen; ++i){
    Rprintf("DEBUG: 12-%d\n", i);
    SEXP thisName = PROTECT(allocVector(STRSXP, 1));
    Rprintf("DEBUG: 13-%d\n", i);
    SET_STRING_ELT(thisName, 0, mkChar(listnames.at(i)));
    Rprintf("DEBUG: 14-%d\n", i);
    SET_VECTOR_ELT(listentrynames, i, thisName);
    Rprintf("DEBUG: 15-%d\n", i);
  }

  SEXP result = PROTECT(allocVector(VECSXP, listlen));
  // 0: tweak / g
  SEXP g_r = PROTECT(allocVector(REALSXP, model.dim()));
  fillPtr(REAL(g_r), model.getParams());
  SET_VECTOR_ELT(result, 0, g_r);
  // 1: Convergence
  SEXP conv_r = PROTECT(allocVector(REALSXP, conv_vec.size()));
  for( std::size_t i = 0; i < conv_vec.size(); ++i){
    Rprintf("from R: conv_r[%d] = %f\n", i, conv_vec(i));
  }
  fillPtr(REAL(conv_r), conv_vec);
  SET_VECTOR_ELT(result, 1, conv_r);
  // 2: lambda:
  SET_VECTOR_ELT(result, 2, _lambda); // simply copy back as the SEXP itself
  setAttrib(result, R_NamesSymbol, listentrynames);
  UNPROTECT(listlen + 4);
  return result;
}
SEXP dtd_evaluate_model_goertler(SEXP model_) {
  auto model = make_model(model_);
  SEXP res = PROTECT(allocVector(REALSXP, 1));
  REAL(res)[0] = model.evaluate();
  UNPROTECT(1);
  return res;
}
extern "C" {
  static const R_CallMethodDef callMethods[] = {
                                                { "_dtd_solve_fista_goertler", (DL_FUNC)&dtd_solve_fista_goertler, 3},
                                                { "_dtd_evaluate_model_goertler", (DL_FUNC)&dtd_evaluate_model_goertler, 1},
                                                {NULL, NULL, 0}
  };
  void R_init_DTD(DllInfo* info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, TRUE);
  }
}
