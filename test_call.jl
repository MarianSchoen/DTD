include("../dtd_julia/main.jl")
using Cxx, Libdl
using StatsBase, Statistics
using Test

# this is meant to be run from within the build directory
addHeaderDir("./src/")
addHeaderDir("/usr/include/eigen3/")
Libdl.dlopen(pwd() * "/build/libdtd_jl.so", Libdl.RTLD_GLOBAL)
cxxinclude(pwd() * "/src/interface_jl.hpp")

const ngenes=100
const ncells=5
const nsamples=10

x = rand(ngenes, ncells)
g = rand(ngenes)
gp = rand(ngenes)
y = rand(ngenes, nsamples)
c = rand(ncells, nsamples)

# define in the cpp space
ptr=pointer(x)
rows=size(x, 1)
cols=size(x, 2)
matcols=rows*cols
cxx_x = icxx"ptrToVec($ptr, $matcols);"

ptr=pointer(y)
rows=size(y, 1)
cols=size(y, 2)
matcols=rows*cols
cxx_y = icxx"ptrToVec($ptr, $matcols);"

ptr=pointer(c)
rows=size(c, 1)
cols=size(c, 2)
matcols=rows*cols
cxx_c = icxx"ptrToVec($ptr, $matcols);"

ptrgp=pointer(gp)
rows = size(gp,1)
cxx_gp = icxx"ptrToVec($ptrgp, $rows);"

ptrg=pointer(g)
rows = size(g,1)
cxx_g = icxx"ptrToVec($ptrg, $rows);"

# TODO: test that matrix conserves rows & cols...

# XTGX_inv:
res=icxx"xtgxinv($cxx_x, $cxx_g, $ngenes, $ncells);"

res_v = collect(res)
xtgxi_cxx = reshape(collect(res), (ncells, ncells));

@test isapprox(xtgxi_cxx, xtgx_inv(x, g))

# test stats function cov, var
res=icxx"cov($cxx_g, $cxx_gp);" # returns a number
@test isapprox(res, cov(g, gp))
res=icxx"var($cxx_g);" # returns a number
@test isapprox(res, var(g))

# gradient:
