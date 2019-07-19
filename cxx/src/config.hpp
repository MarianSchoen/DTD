#pragma once
#include "Eigen/Core"
// TODO: think if single precision is enough? probably.. or first solve on floats, then iterate on double?
typedef double ftype;

// typedef Eigen::Array<ftype, Eigen::Dynamic, 1> vec;
// TODO: for the reference matrix, a sparse matrix type may also be interesting?? Does it pay off?
// typedef Eigen::Array<ftype, Eigen::Dynamic, 2> mat;
using vec = Eigen::VectorXd;
using mat = Eigen::MatrixXd;
