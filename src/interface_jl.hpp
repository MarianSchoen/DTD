#include <iostream> //TODO remove

// this file is the julia interface.



std::vector<double> ptrToVec(double const * ptr, unsigned int len);
std::vector<double> xtgxinv(std::vector<double> const & x, std::vector<double> const & g, std::size_t ngenes, std::size_t ncells);
