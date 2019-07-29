#include "dtd.hpp"
#include "Eigen/Core"
#include <R.h>
#include <Rdefines.h>
#include <array>
#include <string>

// copied from https://github.com/rehbergT/dgemm/blob/master/dgemmR/src/wrapper.cpp
SEXP getElementFromRList(SEXP list_r, const char* name) {
    SEXP names = getAttrib(list_r, R_NamesSymbol);
    for (uint32_t i = 0; i < (uint32_t)length(list_r); i++) {
        if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
            return VECTOR_ELT(list_r, i);
        }
    }
    Rprintf("DTD interface: entry %s not found in list!\n", name);
    return R_NilValue;
}
MatrixXd getMatrixFromR(SEXP mat) {
  double * dmat = REAL(mat);
  int m = INTEGER(GET_DIM(mat))[0];
  int n = INTEGER(GET_DIM(mat))[1];
  return MatrixXd(Eigen::Map<MatrixXd>(dmat, m, n)); // Actually DO make a copy, here (to align memory and be safe)
}
VectorXd getVectorFromR(SEXP v) {
  double* dv = REAL(v);
  int m = XLENGTH(v);
  return VectorXd(Eigen::Map<VectorXd>(dv, m));
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

SEXP dtd_solve_fista_goertler(SEXP model_, SEXP _lambda, SEXP _maxiter, SEXP _saveHistory){
  double lambda = REAL(_lambda)[0];
  int maxiter = INTEGER(_maxiter)[0];
  bool saveHistory = LOGICAL(_saveHistory)[0];
  auto model = make_model(model_);

  dtd::solvers::FistaSolver<dtd::models::GoertlerModel> solver;

  maxiter = std::max(2, maxiter); // "no" error handling

  VectorXd conv_vec(maxiter-1);
  MatrixXd history;
  // first elements are just the status before the iteration:
  conv_vec(0) = model.evaluate();
  if( saveHistory ) {
    history.resize(maxiter-1, model.dim());
    history.row(0) = model.getParams();
  }
  // not the true "iter", but the actual iteration count
  // 0 is the initial value (iter - 1)
  int iter = 1;
  std::function<void(dtd::models::GoertlerModel const & m)> record_solve =
    [&conv_vec,&history,&iter,saveHistory](dtd::models::GoertlerModel const & m) {
      conv_vec(iter) = m.evaluate();
      if( saveHistory ) {
        // this should never happen, but to be sure and not segfault:
        // (and since we cannot throw...)
        if( iter < history.rows() )
          history.row(iter) = m.getParams();
      }
      iter++;
    };

  solver.solve(model, maxiter, lambda, record_solve);

  std::vector<std::string> listnames = {"Tweak", "Convergence", "Lambda"};
  if( saveHistory )
    listnames.push_back("History");
  const std::size_t listlen = listnames.size();
  SEXP listentrynames = PROTECT(allocVector(VECSXP, listlen));
  for( auto i = 0u; i < listlen; ++i){
    SEXP thisName = PROTECT(allocVector(STRSXP, 1));
    SET_STRING_ELT(thisName, 0, mkChar(listnames.at(i).c_str()));
    SET_VECTOR_ELT(listentrynames, i, thisName);
  }

  SEXP result = PROTECT(allocVector(VECSXP, listlen));
  // 0: tweak / g
  SEXP g_r = PROTECT(allocVector(REALSXP, model.dim()));
  fillPtr(REAL(g_r), model.getParams());
  SET_VECTOR_ELT(result, 0, g_r);
  // 1: Convergence
  SEXP conv_r = PROTECT(allocVector(REALSXP, conv_vec.size()));
  fillPtr(REAL(conv_r), conv_vec);
  SET_VECTOR_ELT(result, 1, conv_r);
  // 2: lambda:
  // may also simply copy back as the SEXP itself,
  // but deep-copying makes the UNPROTECT cleaner...
  SEXP _newLambda = PROTECT(allocVector(REALSXP, 1));
  REAL(_newLambda)[0] = lambda;
  SET_VECTOR_ELT(result, 2, _newLambda);
  // 3: History:
  if( saveHistory ) {
    SEXP history_r = PROTECT(allocMatrix(REALSXP, model.dim(), maxiter-1));
    fillPtr(REAL(history_r), history.transpose());
    SET_VECTOR_ELT(result, 3, history_r);
  }

  // set names in list:
  setAttrib(result, R_NamesSymbol, listentrynames);

  // each element in the list has a name and value, hence, 2*listlen,
  // +2 from the list and its vector of names
  UNPROTECT(2*listlen + 2);
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
                                                { "_dtd_solve_fista_goertler", (DL_FUNC)&dtd_solve_fista_goertler, 4},
                                                { "_dtd_evaluate_model_goertler", (DL_FUNC)&dtd_evaluate_model_goertler, 1},
                                                {NULL, NULL, 0}
  };
  void R_init_DTD(DllInfo* info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, TRUE);
  }
}
