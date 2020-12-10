#include "dtd.hpp"
#include "Eigen/Core"
#include <R.h>
#include <Rdefines.h>
#include <array>
#include <string>
#include "nnls.hpp"

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

dtd::models::NormFunctions intToNormFn(int id){
  if( id == 0 ) {
    return dtd::models::NormFunctions::IDENTITY;
  } else if( id == 1 ) {
    return dtd::models::NormFunctions::NORM2;
  } else if( id == 2 ) {
    return dtd::models::NormFunctions::NORM1;
  } else {
    throw std::runtime_error("no such norm function.");
  }
}
dtd::models::ThresholdFunctions intToThreshFn(int id){
  if( id == 0 ) {
    return dtd::models::ThresholdFunctions::SOFTMAX;
  } else {
    throw std::runtime_error("no such threshold function.");
  }
}
dtd::models::SubspaceFunctions intToSubspFn(int id){
  if( id == 0 ) {
    return dtd::models::SubspaceFunctions::POSITIVE;
  } else {
    throw std::runtime_error("no such subspace function.");
  }
}
dtd::models::CoeffEstimation intToCoeffEstim(int id){
  if( id == 0 ) {
    return dtd::models::CoeffEstimation::DIRECT;
  } else if( id == 1 ) {
    return dtd::models::CoeffEstimation::NNLS;
  } else {
    throw std::runtime_error("no such subspace function.");
  }
}

dtd::models::GoertlerModel make_model(SEXP model_) {
  auto x = getMatrixFromR(getElementFromRList(model_, "X"));
  auto y = getMatrixFromR(getElementFromRList(model_, "Y"));
  auto c = getMatrixFromR(getElementFromRList(model_, "C"));
  auto g = getVectorFromR(getElementFromRList(model_, "tweak"));
  int normid = INTEGER(getElementFromRList(model_, "normfnid"))[0];
  int threshfnid = INTEGER(getElementFromRList(model_, "threshfnid"))[0];
  int subspfnid = INTEGER(getElementFromRList(model_, "subspfnid"))[0];
  int estim_c = INTEGER(getElementFromRList(model_, "coeffestim"))[0];
  double inv_prec = REAL(getElementFromRList(model_, "inversion_precision"))[0];

  // crash early (but leave the details to the R space, this is just to prevent segfaults, etc.)
  if( x.size() == 0 || y.size() == 0 || c.size() == 0 || g.size() == 0)
    throw std::runtime_error("make_model: (some of the) input data has zero-length. Cannot build a valid model from it.");
  if( x.rows() !=  g.size() || x.rows() != y.rows() )
    throw std::runtime_error("make_model: x, y and g have a different number of features.");
  if( x.cols() != c.rows() )
    throw std::runtime_error("make_model: x and c have a different number of cells.");
  if( y.cols() != c.cols() )
    throw std::runtime_error("make_model: y and c have a different number of samples.");

  if( inv_prec <= 0 )
    throw std::runtime_error("make_model: inversion_precision is non-positive.");

  return dtd::models::GoertlerModel(x,y,c,g, intToNormFn(normid), intToThreshFn(threshfnid), intToSubspFn(subspfnid), intToCoeffEstim(estim_c), inv_prec);
}

SEXP dtd_solve_fista_goertler(SEXP model_, SEXP _lambda, SEXP _maxiter, SEXP _epsilon, SEXP _navg, SEXP _saveHistory, SEXP _learningrate, SEXP _linesearchspeed, SEXP _cycles, SEXP _restarts, SEXP _haveLearningrate, SEXP _verbose ){
  double lambda = REAL(_lambda)[0];
  int maxiter = INTEGER(_maxiter)[0];
  double epsilon = REAL(_epsilon)[0];
  std::size_t navg = static_cast<std::size_t>(INTEGER(_navg)[0]);
  bool saveHistory = LOGICAL(_saveHistory)[0];
  bool haveLearningrate = LOGICAL(_haveLearningrate)[0];
  double learningrate = 0.0; // will always throw if left uninitialized
  if( haveLearningrate )
    learningrate = REAL(_learningrate)[0];
  double linesearchspeed = REAL(_linesearchspeed)[0];
  std::size_t cycles = static_cast<std::size_t>(INTEGER(_cycles)[0]);
  bool restarts = LOGICAL(_restarts)[0];
  bool verbose = LOGICAL(_verbose)[0];


  try {
    auto model = make_model(model_);
    if( maxiter < 2 )
      throw std::runtime_error("maxiter must be >= 2");

    dtd::solvers::FistaSolver<dtd::models::GoertlerModel> solver;
    if(not haveLearningrate)
      solver.setLearningAuto(model);
    else
      solver.setLearningRate(learningrate);
    solver.setLinesearchSpeed(linesearchspeed);
    solver.setCyclelength(cycles);
    solver.enableRestart(restarts);

    VectorXd conv_vec(maxiter);
    MatrixXd history;
    // first elements are just the status before the iteration:
    if( conv_vec.size() > 0 )
      conv_vec(0) = model.evaluate();
    if( saveHistory ) {
      history.resize(maxiter, model.dim());
      history.row(0) = model.getParams();
    }
    // not the true "iter", but the actual iteration count
    // 0 is the initial value (iter - 1)
    int iter = 1;
    std::function<void(dtd::models::GoertlerModel const &, vec const &)>
        record_solve = [&conv_vec, &history, &iter, saveHistory, verbose,
                        &solver, epsilon](dtd::models::GoertlerModel const &m,
                                          vec const &paramvec) {
          double loss = m.evaluate(paramvec);
          if (iter < conv_vec.size())
            conv_vec(iter) = loss;
          if (saveHistory) {
            // this should never happen, but to be sure and not segfault:
            // (and since we cannot throw...)
            if (iter < history.rows()) {
              history.row(iter) = paramvec;
            }
          }
          if (verbose) {
            Rprintf("**********************************************************"
                    "**********************\n");
            Rprintf("* iteration           %d\n", iter);
            Rprintf("* learning rate:      %5e\n", solver.getLearningRate());
            Rprintf("* loss:               %5f\n", loss);
            Rprintf("* delta y (before N): %5e\n",
                    solver.getDeltaYBeforeNesterov(), epsilon);
            Rprintf("* delta y (after N):  %5e\n",
                    solver.getDeltaYAfterNesterov(), epsilon);
            Rprintf("* nesterov counter:   %d\n", solver.getNesterovCounter());
            Rprintf("* nesterov factor:    %5f\n",
                    dtd::solvers::nesterov_factor(solver.getNesterovCounter()));
          }
          iter++;
        };

    solver.solve(model, maxiter, epsilon, navg, lambda, record_solve);

    if( verbose ) {
      Rprintf("********************************************************************************\n");
      if( iter == maxiter )
        Rprintf("* FISTA did not converge to %e during %d iterations of FISTA\n", epsilon, iter);
      else
        Rprintf("* FISTA successfully converged to %e after %d iterations.\n", epsilon, iter);
      Rprintf("********************************************************************************\n");
    }

    history.conservativeResize(iter, Eigen::NoChange_t());
    conv_vec.conservativeResize(iter);

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
      SEXP history_r = PROTECT(allocMatrix(REALSXP, model.dim(), iter));
      fillPtr(REAL(history_r), history);
      SET_VECTOR_ELT(result, 3, history_r);
    }

    // set names in list:
    setAttrib(result, R_NamesSymbol, listentrynames);

    // each element in the list has a name and value, hence, 2*listlen,
    // +2 from the list and its vector of names
    UNPROTECT(2*listlen + 2);
    return result;
  } catch( std::exception const & exc ){
    error(exc.what());
  }
  return R_NilValue;
}
SEXP dtd_estimate_c(SEXP model_) {
  try {
    auto model = make_model(model_);
    auto c_res = model.estimate_c(model.getParams());
    SEXP res = PROTECT(allocVector(REALSXP, c_res.size()));
    for( int i = 0; i < c_res.size(); ++i ){
      REAL(res)[i] = c_res(i);
    }
    UNPROTECT(1);
    return res;
  } catch( std::exception const & exc ){
    error(exc.what());
  }
  return R_NilValue;
}
SEXP dtd_evaluate_model_goertler(SEXP model_) {
  try {
    auto model = make_model(model_);
    SEXP res = PROTECT(allocVector(REALSXP, 1));
    REAL(res)[0] = model.evaluate();
    UNPROTECT(1);
    return res;
  } catch( std::exception const & exc ){
    error(exc.what());
  }
  return R_NilValue;
}
SEXP dtd_nnls(SEXP mat, SEXP vec, SEXP eps, SEXP maxiter_) {
  try {
    MatrixXd a = getMatrixFromR(mat);
    MatrixXd b = getVectorFromR(vec);
    double epsilon = REAL(eps)[0];
    int maxiter = INTEGER(maxiter_)[0];
    auto c_res = nnls(a, b, epsilon, maxiter);
    SEXP r_res = PROTECT(allocVector(REALSXP, c_res.size()));
    for( int i = 0; i < c_res.size(); ++i) {
      REAL(r_res)[i] = c_res(i);
    }
    UNPROTECT(1);
    return r_res;
  } catch (std::exception const & exc ){
    error(exc.what());
  }
  return R_NilValue;
}
extern "C" {
  static const R_CallMethodDef callMethods[] = {
                                                { "_dtd_solve_fista_goertler", (DL_FUNC)&dtd_solve_fista_goertler, 12},
                                                { "_dtd_evaluate_model_goertler", (DL_FUNC)&dtd_evaluate_model_goertler, 1},
                                                { "_dtd_estimate_c", (DL_FUNC)&dtd_estimate_c, 1},
                                                { "_nnls", (DL_FUNC)&dtd_nnls, 3},
                                                {NULL, NULL, 0}
  };
  void R_init_DTD(DllInfo* info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, TRUE);
  }
}
