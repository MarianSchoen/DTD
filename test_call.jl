include("../dtd_julia/main.jl")
using Cxx, Libdl
# this is meant to be run from within the build directory
addHeaderDir("../src/")
addHeaderDir("/usr/include/eigen3/")
Libdl.dlopen(pwd() * "/libdtd.so", Libdl.RTLD_GLOBAL)
cxxinclude("../src/fista.hpp")

g = rand(100)
x = rand(100,5)

#s = @cxxnew dtd::FistaSolver(pointer(testdata), length(testdata));
g_cxx = @cxxnew vecFromPtr(pointer(g), length(g));
x_cxx = @cxxnew matFromPtr(pointer(x), size(x, 1), size(x, 2))

xinv_cxx = @cxx dtd::invxtgx(x_cxx, g_cxx)
xinv_jl =  xtgx_inv(x, g)

print(xinv_cxx)
print(xinv_jl)

# using Testing
# @test isapprox(xinv_cxx, xinv_jl)

# solver = @cxxnew dtd::FistaSolver(g_cxx);
# @cxx solver->iterate();
