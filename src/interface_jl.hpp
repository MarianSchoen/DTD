#include <iostream> //TODO remove
#include "models.hpp"

// this file is the julia interface.



std::vector<double> ptrToVec(double const * ptr, unsigned int len);
// below is only interfaced for testing. no real usecase for these functions otherwise
std::vector<double> xtgxinv(std::vector<double> const & x, std::vector<double> const & g, std::size_t ngenes, std::size_t ncells);
double var(std::vector<double> const & x);
double cov(std::vector<double> const & x, std::vector<double> const & y);
dtd::models::GoertlerModel makeGoertlerModel(std::vector<double> const & xv, std::vector<double> const & yv, std::vector<double> const & cv, std::size_t ngenes, std::size_t ncells, std::size_t nsamples);
double evalGoertlerModel(dtd::models::GoertlerModel const * model, std::vector<double> const & g);
std::vector<double> gradGoertlerModel(dtd::models::GoertlerModel const * model, std::vector<double> const & g);
