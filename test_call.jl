include("../dtd_julia/main.jl")
using Cxx, Libdl
using StatsBase, Statistics
using Test

# this is meant to be run from within the build directory
addHeaderDir("../src/")
addHeaderDir("/usr/include/eigen3/")
Libdl.dlopen(pwd() * "/libdtd.so", Libdl.RTLD_GLOBAL)
cxxinclude("../src/fista.hpp")

const ngenes=100
const ncells=5
g = rand(ngenes)
h = rand(ngenes)
x = rand(ngenes,ncells)

#s = @cxxnew dtd::FistaSolver(pointer(testdata), length(testdata));
g_cxx = @cxxnew vecFromPtr(pointer(g), length(g));
h_cxx = @cxxnew vecFromPtr(pointer(h), length(h));
x_cxx = @cxxnew matFromPtr(pointer(x), size(x, 1), size(x, 2))

var_cxx = @cxx dtd::stat::var(g_cxx)
var_jl = var(g)
@test isapprox(var_cxx, var_jl)

cov_cxx = @cxx dtd::stat::cov(g_cxx, h_cxx)
cov_jl = cov(g, h)
@test isapprox(cov_cxx, cov_jl)

# this works but is awkward. why does returning a vector not work??
xi = Matrix{Float64}(undef, ncells, ncells)
xi_cxx = @cxxnew matFromPtr(pointer(xi), size(xi, 1), size(xi, 2))
icxx" dtd::invxtgx($x_cxx, $g_cxx, $xi_cxx); "
xi_cxx = reshape(collect(icxx" vecFromMat($xi_cxx); "), ncells, ncells)

@test isapprox(xi_cxx, xtgx_inv(x,g))





# @test isapprox(cov_cxx, cov_jl)

# using Testing
# @test isapprox(xinv_cxx, xinv_jl)

# solver = @cxxnew dtd::FistaSolver(g_cxx);
# @cxx solver->iterate();
